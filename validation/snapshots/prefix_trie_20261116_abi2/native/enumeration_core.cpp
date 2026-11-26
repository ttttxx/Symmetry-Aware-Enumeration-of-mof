// Optional C++17 core. C ABI keeps the Python wrapper dependency-free.
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <exception>
#include <functional>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {
using Clock = std::chrono::steady_clock;
using I = std::int32_t;
using U = std::uint64_t;
thread_local std::string error;
struct Record { std::string key; I stabilizer; };
struct PrefixStats {
    U queries=0,hits=0,early_hits=0,closed=0,created=0,live=0,peak=0,budget_skips=0,snapshots=0;
};
struct PrefixNode {
    U *live;
    bool closed=false;
    std::vector<std::pair<I,std::unique_ptr<PrefixNode>>> children;
    explicit PrefixNode(U *counter):live(counter) { ++*live; }
    ~PrefixNode() { --*live; }
    PrefixNode *find(I value) const {
        for (const auto &child:children) if (child.first==value) return child.second.get();
        return nullptr;
    }
};
struct PrefixForest {
    std::vector<std::vector<I>> partition, vertex_domains;
    std::unique_ptr<PrefixNode> root;
};
struct Context {
    I n;
    std::vector<std::vector<I>> inverses;
    std::unordered_set<std::string> seen;
    std::vector<Record> records;
    bool prefix_enabled=false;
    U prefix_max_nodes=50000;
    PrefixStats prefix_stats;
    // Forests must be destroyed before the counters referenced by their nodes.
    std::vector<std::unique_ptr<PrefixForest>> prefix_forests;
    explicit Context(I sites):n(sites) {}
    std::unique_ptr<PrefixNode> new_prefix_node() {
        if (prefix_stats.live>=prefix_max_nodes) { ++prefix_stats.budget_skips; return nullptr; }
        auto node=std::make_unique<PrefixNode>(&prefix_stats.live);
        ++prefix_stats.created;
        prefix_stats.peak=std::max(prefix_stats.peak,prefix_stats.live);
        return node;
    }
    PrefixNode *prefix_child(PrefixNode *parent,I value) {
        if (auto child=parent->find(value)) return child;
        auto child=new_prefix_node();
        auto result=child.get();
        if (child) parent->children.emplace_back(value,std::move(child));
        return result;
    }
};
struct Task {
    Context &ctx;
    I m, words;
    const I *nodes, *offsets, *domains, *domain_offsets;
    const I *event_offsets, *event_kind, *event_other;
    const U *event_masks, *decided;
    std::int64_t limit;
    bool timed;
    Clock::time_point deadline;
    std::vector<I> coloring, assignment;
    U leaves=0, duplicates=0, pruned=0, candidates=0;
    bool expired=false, capped=false;
    U visits=0;
    std::vector<I> progress{},frontier_values{},frontier_progress{},node_orbit{},expressions{};
    std::vector<std::vector<I>> vertex_domains{};
    std::vector<std::vector<const PrefixForest *>> prefix_queries{};
    I frontier_depth=-1;

    bool stopped(I depth=-1) {
        bool stop=false;
        if (limit >= 0 && ctx.records.size() >= static_cast<U>(limit)) {
            capped=true; stop=true;
        }
        if (timed && Clock::now() >= deadline) { expired=true; stop=true; }
        if (stop && ctx.prefix_enabled && depth>=0 && frontier_depth<0) {
            frontier_depth=depth;
            frontier_values.assign(assignment.begin(),assignment.begin()+depth);
            frontier_progress=progress;
        }
        return stop;
    }
    void prepare_prefix_queries() {
        if (!ctx.prefix_enabled) return;
        progress.assign(m,0);
        node_orbit.assign(ctx.n,0);
        expressions.assign(ctx.n,0);
        prefix_queries.resize(m+1);
        vertex_domains.resize(ctx.n);
        for (I i=0;i<m;++i) {
            std::vector<I> choices(domains+domain_offsets[i],domains+domain_offsets[i+1]);
            std::sort(choices.begin(),choices.end());
            choices.erase(std::unique(choices.begin(),choices.end()),choices.end());
            for (I j=offsets[i];j<offsets[i+1];++j) {
                I node=nodes[j];
                node_orbit[node]=i;
                vertex_domains[node]=choices;
            }
        }
        for (const auto &forest:ctx.prefix_forests) {
            bool subset=true;
            for (I node=0;node<ctx.n && subset;++node) {
                const auto &source=forest->vertex_domains[node],&current=vertex_domains[node];
                subset=std::includes(source.begin(),source.end(),current.begin(),current.end());
            }
            if (!subset) continue;
            I equality_depth=0;
            bool impossible=false;
            for (const auto &orbit:forest->partition) {
                std::vector<I> terms;
                for (I node:orbit) {
                    I variable=node_orbit[node];
                    I term=domain_offsets[variable+1]-domain_offsets[variable]==1 ?
                        domains[domain_offsets[variable]] : -variable-1;
                    if (std::find(terms.begin(),terms.end(),term)==terms.end()) terms.push_back(term);
                }
                if (terms.size()<=1) continue;
                I constants=0;
                for (I term:terms) {
                    if (term>=0) ++constants;
                    else equality_depth=std::max(equality_depth,-term);
                }
                if (constants>1) { impossible=true;break; }
            }
            if (impossible) continue;
            struct Pending {const PrefixNode *node;I level,required;};
            std::vector<Pending> stack{{forest->root.get(),0,equality_depth}};
            std::vector<bool> triggers(m+1,false);
            while (!stack.empty()) {
                auto frame=stack.back();stack.pop_back();
                if (frame.node->closed) { triggers[frame.required]=true;continue; }
                if (frame.node->children.empty()) continue;
                I variable=node_orbit[forest->partition[frame.level][0]];
                I known=domain_offsets[variable+1]-domain_offsets[variable]==1 ? 0 : variable+1;
                for (const auto &child:frame.node->children)
                    stack.push_back({child.second.get(),frame.level+1,std::max(frame.required,known)});
            }
            for (I depth=0;depth<=m;++depth)
                if (triggers[depth]) prefix_queries[depth].push_back(forest.get());
        }
    }
    bool prefix_covered(I depth) {
        if (!ctx.prefix_enabled || prefix_queries[depth].empty()) return false;
        ++ctx.prefix_stats.queries;
        for (I node=0;node<ctx.n;++node) {
            I variable=node_orbit[node];
            expressions[node]=coloring[node]>=0 ? coloring[node] :
                (domain_offsets[variable+1]-domain_offsets[variable]==1 ?
                 domains[domain_offsets[variable]] : -variable-1);
        }
        for (const auto *forest:prefix_queries[depth]) {
            bool fixed=true;
            for (const auto &orbit:forest->partition) {
                I first=expressions[orbit[0]];
                for (I node:orbit) if (expressions[node]!=first) { fixed=false; break; }
                if (!fixed) break;
            }
            if (!fixed) continue;
            const PrefixNode *current=forest->root.get();
            for (const auto &orbit:forest->partition) {
                if (current->closed) break;
                I value=expressions[orbit[0]];
                if (value<0) { current=nullptr; break; }
                current=current->find(value);
                if (!current) break;
            }
            if (current && current->closed) {
                ++ctx.prefix_stats.hits;
                if (depth<m) ++ctx.prefix_stats.early_hits;
                return true;
            }
        }
        return false;
    }
    void capture_frontier() {
        if (!ctx.prefix_enabled || frontier_depth<0) return;
        ++ctx.prefix_stats.snapshots;
        auto forest=std::make_unique<PrefixForest>();
        forest->root=ctx.new_prefix_node();
        if (!forest->root) return;
        forest->vertex_domains=vertex_domains;
        for (I i=0;i<m;++i) forest->partition.emplace_back(nodes+offsets[i],nodes+offsets[i+1]);
        PrefixNode *current=forest->root.get();
        U closed=0;
        for (I level=0;level<=frontier_depth && level<m;++level) {
            for (I choice=0;choice<frontier_progress[level];++choice) {
                I value=domains[domain_offsets[level]+choice];
                PrefixNode *child=ctx.prefix_child(current,value);
                if (child && !child->closed) {
                    child->children.clear();child->closed=true;++closed;
                }
            }
            if (level>=frontier_depth || current->closed) break;
            current=ctx.prefix_child(current,frontier_values[level]);
            if (!current || current->closed) break;
        }
        if (closed) {
            ctx.prefix_stats.closed+=closed;
            ctx.prefix_forests.push_back(std::move(forest));
        }
    }
    bool excluded(I depth, std::vector<U> &mask) {
        if (depth) {
            I index=depth-1;
            for (I e=event_offsets[index]; e<event_offsets[index+1]; ++e) {
                I other=event_kind[e] == 0 ? event_other[e] : assignment[event_other[e]];
                if (assignment[index] != other)
                    for (I w=0; w<words; ++w) mask[w] &= ~event_masks[e*words+w];
            }
        }
        for (I w=0; w<words; ++w)
            if (mask[w] & decided[depth*words+w]) { ++pruned; return true; }
        return false;
    }
    void accept(I stabilizer) {
        ++leaves;
        std::string key(ctx.n, '\0');
        for (I i=0; i<ctx.n; ++i) key[i]=static_cast<char>(coloring[i]);
        if (!ctx.seen.insert(key).second) { ++duplicates; return; }
        ctx.records.push_back({std::move(key), stabilizer});
    }
    bool dfs(I depth, const std::vector<std::pair<I,I>> &active, std::vector<U> mask) {
        ++visits;
        if (ctx.prefix_enabled && depth<m) progress[depth]=0;
        if (stopped(depth)) return false;
        if (excluded(depth, mask) || prefix_covered(depth)) return true;
        if (depth == m) { accept(static_cast<I>(active.size())+1); return true; }
        for (I d=domain_offsets[depth]; d<domain_offsets[depth+1]; ++d) {
            if (stopped(depth)) return false;
            I value=domains[d];
            assignment[depth]=value;
            for (I j=offsets[depth]; j<offsets[depth+1]; ++j) coloring[nodes[j]]=value;
            std::vector<std::pair<I,I>> next;
            next.reserve(active.size());
            bool reject=false;
            for (auto [p, cursor]: active) {
                const auto &inverse=ctx.inverses[p];
                while (cursor < ctx.n) {
                    I current=coloring[cursor], image=coloring[inverse[cursor]];
                    if (current < 0 || image < 0) { next.emplace_back(p,cursor); break; }
                    if (image < current) { reject=true; break; }
                    if (image > current) break;
                    ++cursor;
                }
                if (cursor == ctx.n) next.emplace_back(p,cursor);
                if (reject) break;
            }
            bool completed=reject || dfs(depth+1, next, mask);
            for (I j=offsets[depth]; j<offsets[depth+1]; ++j) coloring[nodes[j]]=-1;
            if (!completed) return false;
            if (ctx.prefix_enabled) progress[depth]=d-domain_offsets[depth]+1;
        }
        return true;
    }
    void counter(const U *initial) {
        std::vector<I> digit(m,0), node_orbit(ctx.n);
        for (I i=0; i<m; ++i) {
            if (domain_offsets[i] == domain_offsets[i+1]) return;
            for (I j=offsets[i]; j<offsets[i+1]; ++j) node_orbit[nodes[j]]=i;
        }
        std::vector<std::vector<std::pair<I,I>>> pairs;
        for (const auto &inverse:ctx.inverses) {
            std::vector<std::pair<I,I>> compared;
            for (I site=0; site<ctx.n; ++site) {
                I a=node_orbit[site], b=node_orbit[inverse[site]];
                if (a==b) continue;
                if (domain_offsets[a+1]-domain_offsets[a] == 1 &&
                    domain_offsets[b+1]-domain_offsets[b] == 1 &&
                    domains[domain_offsets[a]] == domains[domain_offsets[b]]) continue;
                compared.emplace_back(a,b);
            }
            pairs.push_back(std::move(compared));
        }
        while (!stopped()) {
            ++candidates;
            for (I i=0; i<m; ++i) assignment[i]=domains[domain_offsets[i]+digit[i]];
            std::vector<U> mask(initial, initial+words);
            bool reject=excluded(0,mask);
            for (I depth=1; depth<=m && !reject; ++depth) reject=excluded(depth,mask);
            I stabilizer=1;
            if (!reject) {
                for (const auto &compared:pairs) {
                    bool equal=true;
                    for (auto [a,b]:compared) {
                        if (assignment[b] < assignment[a]) { reject=true; equal=false; break; }
                        if (assignment[b] > assignment[a]) { equal=false; break; }
                    }
                    if (reject) break;
                    if (equal) ++stabilizer;
                }
            }
            if (!reject) {
                for (I i=0; i<ctx.n; ++i) coloring[i]=assignment[node_orbit[i]];
                accept(stabilizer);
            }
            I i=m-1;
            for (; i>=0; --i) {
                if (++digit[i] < domain_offsets[i+1]-domain_offsets[i]) break;
                digit[i]=0;
            }
            if (i<0) return;
        }
    }
};
}

extern "C" {
const char *mof_core_error() { return error.c_str(); }
I mof_core_version() { return 2; }
void *mof_core_create(I n, I group_size, const I *inverse_data) {
    try {
        error.clear();
        if (n<=0 || n>512 || group_size<=0) throw std::invalid_argument("Invalid native dimensions");
        auto ctx=new Context(n);
        try {
            I identities=0;
            for (I p=0; p<group_size; ++p) {
                std::vector<I> inverse(inverse_data+p*n,inverse_data+(p+1)*n);
                std::vector<bool> seen(n,false);
                bool identity=true;
                for (I i=0; i<n; ++i) {
                    if (inverse[i]<0 || inverse[i]>=n || seen[inverse[i]])
                        throw std::invalid_argument("Invalid inverse permutation");
                    seen[inverse[i]]=true;
                    identity &= inverse[i]==i;
                }
                if (identity) ++identities; else ctx->inverses.push_back(std::move(inverse));
            }
            if (identities!=1) throw std::invalid_argument("Expected exactly one identity");
        } catch (...) { delete ctx; throw; }
        return ctx;
    } catch (const std::exception &e) { error=e.what(); return nullptr; }
}
void mof_core_destroy(void *context) { delete static_cast<Context *>(context); }
I mof_core_prefix_mode(void *context,I enabled,U max_nodes) {
    auto &ctx=*static_cast<Context *>(context);
    ctx.prefix_enabled=enabled!=0;ctx.prefix_max_nodes=max_nodes;
    return 0;
}
void mof_core_prefix_stats(void *context,U *out) {
    const auto &ctx=*static_cast<Context *>(context);
    const auto &s=ctx.prefix_stats;
    U values[]={s.queries,s.hits,s.early_hits,s.closed,s.created,s.live,s.peak,
                s.budget_skips,s.snapshots,static_cast<U>(ctx.prefix_forests.size())};
    std::copy(values,values+10,out);
}
std::uint64_t mof_core_count(void *context) { return static_cast<Context *>(context)->records.size(); }
I mof_core_run(void *context, I m, const I *nodes, const I *offsets,
               const I *domains, const I *domain_offsets,
               I words, const U *initial, const U *decided,
               const I *event_offsets, const I *event_kind, const I *event_other,
               const U *event_masks, std::int64_t limit, double remaining,
               I kernel, U *stats) {
    try {
        error.clear();
        auto &ctx=*static_cast<Context *>(context);
        if (m<=0 || m>ctx.n || words<=0 || (kernel!=0 && kernel!=1) || (kernel==1 && ctx.prefix_enabled))
            throw std::invalid_argument("Invalid native task");
        Task task{ctx,m,words,nodes,offsets,domains,domain_offsets,event_offsets,
                  event_kind,event_other,event_masks,decided,limit,remaining>=0,
                  Clock::now()+std::chrono::duration_cast<Clock::duration>(std::chrono::duration<double>(remaining>=0?remaining:0)),
                  std::vector<I>(ctx.n,-1),std::vector<I>(m,-1)};
        if (kernel==0) {
            std::vector<std::pair<I,I>> active;
            for (I p=0; p<static_cast<I>(ctx.inverses.size()); ++p) active.emplace_back(p,0);
            task.prepare_prefix_queries();
            bool completed=task.dfs(0,active,std::vector<U>(initial,initial+words));
            if (!completed) task.capture_frontier();
        } else task.counter(initial);
        stats[0]=task.leaves; stats[1]=task.duplicates; stats[2]=task.pruned;
        stats[3]=task.candidates; stats[4]=task.expired; stats[5]=task.capped;
        stats[6]=task.visits;
        return 0;
    } catch (const std::exception &e) { error=e.what(); return -1; }
}
I mof_core_copy(void *context, U start, U count, unsigned char *keys, I *stabilizers) {
    try {
        auto &ctx=*static_cast<Context *>(context);
        if (start>ctx.records.size() || count>ctx.records.size()-start)
            throw std::invalid_argument("Invalid result slice");
        for (U i=0; i<count; ++i) {
            const auto &record=ctx.records[start+i];
            std::memcpy(keys+i*ctx.n,record.key.data(),ctx.n);
            stabilizers[i]=record.stabilizer;
        }
        return 0;
    } catch (const std::exception &e) { error=e.what(); return -1; }
}
}

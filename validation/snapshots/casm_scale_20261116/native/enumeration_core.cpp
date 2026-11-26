// Optional C++17 core. C ABI keeps the Python wrapper dependency-free.
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <exception>
#include <functional>
#include <limits>
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
struct Context {
    I n;
    std::vector<std::vector<I>> inverses;
    std::unordered_set<std::string> seen;
    std::vector<Record> records;
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

    bool stopped() {
        if (limit >= 0 && ctx.records.size() >= static_cast<U>(limit)) {
            capped=true; return true;
        }
        if (timed && Clock::now() >= deadline) { expired=true; return true; }
        return false;
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
    void dfs(I depth, const std::vector<std::pair<I,I>> &active, std::vector<U> mask) {
        if (stopped() || excluded(depth, mask)) return;
        if (depth == m) { accept(static_cast<I>(active.size())+1); return; }
        for (I d=domain_offsets[depth]; d<domain_offsets[depth+1]; ++d) {
            if (stopped()) return;
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
            if (!reject) dfs(depth+1, next, mask);
            for (I j=offsets[depth]; j<offsets[depth+1]; ++j) coloring[nodes[j]]=-1;
        }
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
I mof_core_version() { return 1; }
void *mof_core_create(I n, I group_size, const I *inverse_data) {
    try {
        error.clear();
        if (n<=0 || n>512 || group_size<=0) throw std::invalid_argument("Invalid native dimensions");
        auto ctx=new Context{n,{},{},{}};
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
        if (m<=0 || m>ctx.n || words<=0 || (kernel!=0 && kernel!=1))
            throw std::invalid_argument("Invalid native task");
        Task task{ctx,m,words,nodes,offsets,domains,domain_offsets,event_offsets,
                  event_kind,event_other,event_masks,decided,limit,remaining>=0,
                  Clock::now()+std::chrono::duration_cast<Clock::duration>(std::chrono::duration<double>(remaining>=0?remaining:0)),
                  std::vector<I>(ctx.n,-1),std::vector<I>(m,-1)};
        if (kernel==0) {
            std::vector<std::pair<I,I>> active;
            for (I p=0; p<static_cast<I>(ctx.inverses.size()); ++p) active.emplace_back(p,0);
            task.dfs(0,active,std::vector<U>(initial,initial+words));
        } else task.counter(initial);
        stats[0]=task.leaves; stats[1]=task.duplicates; stats[2]=task.pruned;
        stats[3]=task.candidates; stats[4]=task.expired; stats[5]=task.capped;
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

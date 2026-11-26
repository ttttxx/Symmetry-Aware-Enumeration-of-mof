// Standalone structural/oracle test; no extra production ABI for test insertion.
#include "../native/enumeration_core.cpp"
#include <cassert>
#include <random>
#include <unordered_map>

int main() {
    for (I length:{1,2,16,512}) {
        std::mt19937 random(47+length);
        std::vector<std::string> inputs;
        for (I k=0;k<2000;++k) {
            std::string key(length,'\0');
            if (k<1000) { // Long common paths plus every unsigned byte code.
                key[length-1]=static_cast<char>(k%256);
                if (length>1) key[length-2]=static_cast<char>(k/256);
            } else for (auto &value:key) value=static_cast<char>(random()%256);
            inputs.push_back(key);
        }
        for (I order=0;order<6;++order) {
            if (order==0) std::sort(inputs.begin(),inputs.end());
            else if (order==1) std::reverse(inputs.begin(),inputs.end());
            else std::shuffle(inputs.begin(),inputs.end(),random);
            CompactTrie trie;
            std::vector<Record> records;
            std::unordered_map<std::string,U> oracle;
            for (const auto &key:inputs) {
                bool expected=oracle.emplace(key,records.size()).second;
                assert(trie.insert(key,records)==expected);
                if (expected) records.push_back({key,1});
                assert(!trie.insert(key,records));
                assert(trie.lookup(reinterpret_cast<const unsigned char *>(key.data()),records)==
                       static_cast<std::int64_t>(oracle.at(key)));
            }
            for (const auto &item:oracle)
                assert(trie.lookup(reinterpret_cast<const unsigned char *>(item.first.data()),records)==
                       static_cast<std::int64_t>(item.second));
            for (I trial=0;trial<1000;++trial) {
                std::string key(length,'\0');
                for (auto &value:key) value=static_cast<char>(random()%256);
                auto found=oracle.find(key);
                auto expected=found==oracle.end() ? -1 : static_cast<std::int64_t>(found->second);
                assert(trie.lookup(reinterpret_cast<const unsigned char *>(key.data()),records)==expected);
            }
            std::vector<bool> visited(trie.nodes.size(),false);
            std::vector<CompactTrie::Index> stack{0};
            while (!stack.empty()) {
                auto parent=stack.back();stack.pop_back();
                assert(!visited[parent]);visited[parent]=true;
                I previous=-1,children=0;
                for (auto child=trie.nodes[parent].child;child!=CompactTrie::nil;child=trie.nodes[child].next) {
                    const auto &node=trie.nodes[child];
                    assert(node.sample<records.size());
                    assert(node.depth>trie.nodes[parent].depth && node.depth<=static_cast<U>(length));
                    I label=static_cast<unsigned char>(records[node.sample].key[trie.nodes[parent].depth]);
                    assert(label>previous);previous=label;
                    ++children;stack.push_back(child);
                }
                if (parent && children) assert(children>=2);
                if (parent && !children) assert(trie.nodes[parent].depth==static_cast<U>(length));
            }
            assert(std::all_of(visited.begin(),visited.end(),[](bool value){return value;}));
            assert(trie.nodes.size()<=2*records.size()+1);
        }
    }
}

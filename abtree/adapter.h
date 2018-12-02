/**
 * Implementation of a lock-free relaxed (a,b)-tree using LLX/SCX.
 * Trevor Brown, 2018.
 */

#ifndef DS_ADAPTER_H
#define DS_ADAPTER_H

#include <iostream>
#include "errors.h"
#include "random_fnv1a.h"
#include "brown_ext_abtree_lf_impl.h"
#ifdef USE_TREE_STATS
#   include "tree_stats.h"
#endif

#if !defined FAT_NODE_DEGREE
    #define FAT_NODE_DEGREE 11
#endif

#define RECORD_MANAGER_T record_manager<Reclaim, Alloc, Pool, abtree_ns::Node<FAT_NODE_DEGREE, K>>
#define DATA_STRUCTURE_T abtree_ns::abtree<FAT_NODE_DEGREE, K, std::less<K>, RECORD_MANAGER_T>

template <typename K, typename V, class Reclaim = reclaimer_debra<K>, class Alloc = allocator_new<K>, class Pool = pool_none<K>>
class ds_adapter {
private:
    DATA_STRUCTURE_T * const ds;

public:
    ds_adapter(const int NUM_THREADS,
               const K& KEY_ANY,
               const K& unused1,
               const V& unused2,
               RandomFNV1A * const unused3)
    : ds(new DATA_STRUCTURE_T(NUM_THREADS, KEY_ANY))
    {
        if (!std::is_same<V, void *>::value) {
            setbench_error("Value type V used with brown_ext_abtree_lf is not void *. This data structure has hard writes value type void *.");
        }
        if (NUM_THREADS > MAX_THREADS_POW2) {
            setbench_error("NUM_THREADS exceeds MAX_THREADS_POW2");
        }
    }
    ~ds_adapter() {
        delete ds;
    }
    
    void * getNoValue() {
        return ds->NO_VALUE;
    }
    
    void initThread(const int tid) {
        ds->initThread(tid);
    }
    void deinitThread(const int tid) {
        ds->deinitThread(tid);
    }

    bool contains(const int tid, const K& key) {
        return ds->contains(tid, key);
    }
    void * insert(const int tid, const K& key, void * const val) {
        return ds->insert(tid, key, val);
    }
    void * insertIfAbsent(const int tid, const K& key, void * const val) {
        return ds->insertIfAbsent(tid, key, val);
    }
    void * erase(const int tid, const K& key) {
        return ds->erase(tid, key).first;
    }
    void * find(const int tid, const K& key) {
        return ds->find(tid, key).first;
    }
    int rangeQuery(const int tid, const K& lo, const K& hi, K * const resultKeys, void ** const resultValues) {
        return ds->rangeQuery(tid, lo, hi, resultKeys, resultValues);
    }
    void printSummary() {
        ds->debugGetRecMgr()->printStatus();
    }
    bool validateStructure() {
        return true;
    }
    void printObjectSizes() {
        std::cout<<"size_node="<<(sizeof(abtree_ns::Node<FAT_NODE_DEGREE, K>))<<std::endl;
    }
#ifdef USE_TREE_STATS
    class NodeHandler {
    public:
        typedef abtree_ns::Node<FAT_NODE_DEGREE,K> * NodePtrType;

        K minKey;
        K maxKey;
        
        NodeHandler(const K& _minKey, const K& _maxKey) {
            minKey = _minKey;
            maxKey = _maxKey;
        }
        
        class ChildIterator {
        private:
            size_t ix;
            NodePtrType node; // node being iterated over
        public:
            ChildIterator(NodePtrType _node) { node = _node; ix = 0; }
            bool hasNext() { return ix < node->size; }
            NodePtrType next() { return node->ptrs[ix++]; }
        };
        
        static bool isLeaf(NodePtrType node) { return node->leaf; }
        static ChildIterator getChildIterator(NodePtrType node) { return ChildIterator(node); }
        static size_t getNumChildren(NodePtrType node) { return node->size; }
        static size_t getNumKeys(NodePtrType node) { return isLeaf(node) ? node->size : 0; }
        static size_t getSumOfKeys(NodePtrType node) {
            size_t sz = getNumKeys(node);
            size_t result = 0;
            for (size_t i=0;i<sz;++i) {
                result += (size_t) node->keys[i];
            }
            return result;
        }
    };
    TreeStats<NodeHandler> * createTreeStats(const K& _minKey, const K& _maxKey) {
        return new TreeStats<NodeHandler>(new NodeHandler(_minKey, _maxKey), ds->debug_getEntryPoint(), true);
    }
#endif
};

#endif

/*
 * gbst.h
 *
 * GreenBST
 *
 * This is part of the tree library
 *
 * Copyright 2015 Ibrahim Umar (UiT the Arctic University of Norway)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * ---
 *
 * Several vEB related operations are based on the code from:
 * G. S. Brodal, R. Fagerberg, and R. Jacob, “Cache oblivious search trees via binary trees of small height,”
 * in Proceedings of the thirteenth annual ACM-SIAM symposium on Discrete algorithms, ser. SODA ’02, 2002, pp. 39–48.
 *
 */



#ifndef gbst_h
#define gbst_h

#include <stdint.h>
#include "gbstlock.h"
#include "locks.h"

#define atomic_inc(P) __sync_add_and_fetch((P), 1)
#define atomic_dec(P) __sync_add_and_fetch((P), -1)

//To count the number of waiting thread
//#define __WAIT_COUNT

#ifdef __WAIT_COUNT

unsigned __wait;

static inline int pthread_spin_lock_w(gbst_lock_t *lock)
{
	if (pthread_spin_trylock(lock)) {
		atomic_inc(&__wait);
		pthread_spin_lock(lock);
	}
	return 0;
}

#else

#define pthread_spin_lock_w(x)  pthread_spin_lock(x)

#endif /* __WAIT_COUNT */


#ifdef DEBUG
#define DEBUG_PRINT(...) { fprintf(stderr, "\n>>DEBUG: "); fprintf(stderr, __VA_ARGS__); }
#else
#define DEBUG_PRINT(...) {}
#endif

#define MAXINT  INT_MAX
#define EMPTY   0
#define FULL    1

#define MAX_STACK 50


/* For 32 bit keys (unsigned) */

#define _NODETYPE uint32_t
#define _NODEWIDTH 32
#define _NODESHIFT 2

/* For 64 bit keys (unsigned long) */
/*
#define _NODETYPE uint64_t
#define _NODEWIDTH 64
#define _NODESHIFT 3
*/


struct map { size_t left; size_t right; };

struct node {
	_NODETYPE value;
};


struct GNode {
	unsigned	rev;
	_NODETYPE	high_key;
	struct GNode *sibling;
	unsigned	isleaf;

	gbst_lock_t	lock;           //Tree lock
	unsigned	count_node;     //sub-tree node count
#ifdef __STATS
	unsigned	deleted_node;
#endif
#ifdef __PREALLOCGNODES
	struct node	a[__PREALLOCGNODES];
	uintptr_t 	b[__PREALLOCGNODES];
#else
	struct node	*a;
	uintptr_t 	*b;
#endif
};


struct global {
	int		max_node;
	int		max_depth;
	int		split_thres;

#ifdef __STATS
	int		nodecnt;
	int		maxcnt;

	int		count_ins;
	int		count_del;

	int		failed_ins;
	int		failed_del;

	int		rebalance_done_ins;
	int		rebalance_done_del;
#endif
	gbst_lock_t	lock;

	struct GNode ** root;
};




//--STACK START
struct stack {
	int		num;
	struct GNode *	parent[MAX_STACK];
};

static inline
int push(struct stack *stk, void *parent)
{
	stk->parent[stk->num++] = (struct GNode *)parent;
	return 0;
}
static inline
struct GNode *pop(struct stack *stk)
{
	if (stk->num > 0) {
		stk->num--;
		return stk->parent[stk->num];
	} else {
		return 0;
	}
}
//--STACK END


//Map
extern struct map *_map;

#ifdef __PREALLOCGNODES

#define MAX_POOLSIZE 40000

extern struct map mapcontent[__PREALLOCGNODES];

//Pool
extern struct GNode _pool [MAX_POOLSIZE];

//Counter
extern unsigned _poolCtr;

#endif

static inline
_NODETYPE _val(_NODETYPE x)
{
	return x & ~((_NODETYPE)1 << (_NODEWIDTH-1));
}

static inline
void set_mark(_NODETYPE *value)
{
	*value = *value | (_NODETYPE)1 << (_NODEWIDTH-1);
}

static inline
void unset_mark(_NODETYPE *value)
{
	*value = *value & ~((_NODETYPE)1 << (_NODEWIDTH-1));
}

static inline
int is_marked(_NODETYPE value)
{
	return value >> (_NODEWIDTH-1);
}

static inline
int get_idx(void* p, void* base, size_t nodesize)
{
    return (int)(((uintptr_t)p - (uintptr_t)base)/nodesize);
}

/* Below is left and right WITHOUT branching. Careful with the node size. */

//#define left(p,base) (void*) (((uintptr_t)base +_map[((uintptr_t)p - (uintptr_t)base)>>_NODESHIFT].left) * !((SIZE_MAX + _map[((uintptr_t)p - (uintptr_t)base)>>_NODESHIFT].left) >> (_NODEWIDTH-1)))
//#define right(p,base) (void*) (((uintptr_t)base +_map[((uintptr_t)p - (uintptr_t)base)>>_NODESHIFT].right) * !((SIZE_MAX + _map[((uintptr_t)p - (uintptr_t)base)>>_NODESHIFT].right) >> (_NODEWIDTH-1)))

/* Below is left and right WITH branching. Careful with the node size. */

//#define left(p,base) (void*)(_map[(int)(((uintptr_t)p - (uintptr_t)base)>>_NODESHIFT)].left?((uintptr_t)base + _map[(int)(((uintptr_t)p - (uintptr_t)base)>>_NODESHIFT)].left):0)
//#define right(p,base) (void*)(_map[(int)(((uintptr_t)p - (uintptr_t)base)>>_NODESHIFT)].right?((uintptr_t)base + _map[(int)(((uintptr_t)p - (uintptr_t)base)>>_NODESHIFT)].right):0)


static inline
struct node *left(void *p, void *base)
{
	int idx = (int)(((size_t)p - (size_t)base) >> _NODESHIFT);

	//fprintf(stderr, "going left from index %d, to %ld (%d)\n", idx, _map[idx].left + base, _map[idx].left );

	if (_map[idx].left != 0)
		return (struct node *) ((uintptr_t) base + (uintptr_t) _map[idx].left);
	else
		return 0;
}

static inline
struct node *right(void *p, void *base)
{
	int idx = (int)(((size_t)p - (size_t)base) >> _NODESHIFT);

	//fprintf(stderr, "going right from index %d, to %ld (%d)\n", idx, _map[idx].right + base, _map[idx].right );

	if (_map[idx].right != 0)
		return (struct node *) ((uintptr_t) base + (uintptr_t) _map[idx].right);
	else
		return 0;
}

typedef struct global greenbst_t;

greenbst_t *greenbst_alloc(int t);

int greenbst_free(greenbst_t *map);
int greenbst_insert(greenbst_t *map, _NODETYPE key, void *data);
int greenbst_contains(greenbst_t *map, _NODETYPE key);
void *greenbst_get(greenbst_t *map, _NODETYPE key);
int greenbst_delete(greenbst_t *map, _NODETYPE key);

void initial_add(struct global *universe, int num, int range);

#ifndef __PREALLOCGNODES
void init_threads(int);
#endif

#endif /* gbst_h */

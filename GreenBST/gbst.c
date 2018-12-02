/*
 * gbst.c
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

#define _GNU_SOURCE
#include <sched.h>

#include <unistd.h>

#include <assert.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <time.h>
#include <sys/times.h>
#include <sys/time.h>

#include <signal.h>
#include <limits.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <inttypes.h>

#include "gbst.h"

#include "locks.h"
#include "gbstlock.h"
#include <pthread.h>

//---START RECURSIVE BRODAL!!!

typedef struct li {
	int p, ts, bs;
} levelinfo;

void init_height(int top, int bot, levelinfo *myli)
{
	int me = (top + bot + 1) / 2;

	if (top <= bot + 1) return;
	myli[me].p = top;
	myli[me].ts = (1 << (top - me)) - 1;
	myli[me].bs = (1 << (me - bot)) - 1;

	//printf("me: %d, p: %d, ts: %d, bs: %d\n", me, myli[me].p, myli[me].ts, myli[me].bs);

	init_height(top, me, myli);
	init_height(me, bot, myli);
}


void impl_report_al_rec(int ad, int parent, int *val, int *size, int *h, levelinfo *myli, int*depth, int *bf, int *anc)
{
	int left, right;

	if (!(ad < *size && *h > 0)) return;

	//_map[ad].parent = parent;
	//_map[ad].depth = depth;

	*h = *h - 1;
	*depth = *depth + 1;

	left = (*bf = *bf * 2, anc[*h] = (myli[*h].bs * (*bf & myli[*h].ts)) + myli[*h].ts + anc[myli[*h].p]);

	if (*bf > *size)
		_map[ad].left = 0;
	else
		_map[ad].left = left * sizeof(_NODETYPE);

	//printf("No:%d, left:%d, bf:%d, anc[h]:%d \n", ad, left, bf, anc[h]);

	impl_report_al_rec(left, ad, val, size, h, myli, depth, bf, anc);

	*bf = *bf >> 1;



/*
 *
 *      fprintf(stderr,"<%2d>",depth);
 *      for(e=0;e<depth-1;e++) fprintf(stderr," .");
 *      if(val[ ad ]== MAXINT ) fprintf(stderr,"  +");
 *      else fprintf(stderr,"%3d",val[ ad ]);
 *      for(;e<16;e++) fprintf(stderr," -");fprintf(stderr,"<%d>\n",  ad );
 */


	right = (*bf = 2 * *bf + 1, anc[*h] = (myli[*h].bs * (*bf & myli[*h].ts)) + myli[*h].ts + anc[myli[*h].p]);

	if (*bf > *size)
		_map[ad].right = 0;
	else
		_map[ad].right = right * sizeof(_NODETYPE);

	//printf("No:%d, right:%d, bf:%d, anc[h]:%d \n", ad, right, bf, anc[h]);

	impl_report_al_rec(right, ad, val, size, h, myli, depth, bf, anc);

	*h = *h + 1;
	*bf = *bf >> 1;
	*depth = *depth - 1;


	return;
}

void impl_report_al(int *val, int *size, int *h, levelinfo *myli, int*depth, int *bf, int *anc, int max_depth)
{
	*h = max_depth;

	*bf = 1;

	anc[max_depth] = 0;

	impl_report_al_rec(0, 0, val, size, h, myli, depth, bf, anc);
}

void init_tree_map(struct global *universe, int max_node, int max_depth)
{
	levelinfo *myli = 0;

	int depth = 0;

	int size;

	int *anc = 0;

	int bf;

	int h;

#ifdef __PREALLOCGNODES
	_map = &mapcontent[0];
#else
	_map = (struct map *)calloc(max_node, sizeof(struct map));
#endif

	size = max_node;

	if (!myli) myli = (levelinfo *) malloc((max_depth + 2) * sizeof(levelinfo));
	if (!myli) exit(23);

	myli[0].p = 0;
	myli[0].ts = 0;
	myli[0].bs = 0;

	if (!anc) anc = (int*) malloc((max_depth + 2) * sizeof(int));
	if (!anc) exit(23);

	printf("max_depth: %d, %d\n", max_depth, max_node);

	init_height(max_depth, 0, myli);

	//if(val) free(val);
	int *val = (int*) malloc((size + 2) * sizeof(int));
	if (!val) exit(23);

	impl_report_al(val, &size, &h, myli, &depth, &bf, anc, max_depth);

	free(val);
	free(anc);
	free(myli);

	//int d = 0;

	//for (d = 0; d <max_node ; d++)
	//    printf( "<%d>, left %d, right %d\n", d, _map[d].left, _map[d].right);

}

//---END RECURSIVE BRODAL!!!


void draw_helper(struct node *p, int *depth, void *base, int currID)
{
	int i;

	if (p) {
		(*depth)++;
		//if(p->sleft){
		//    struct metadata_struct *mt = ss_lookup(p);
		//    draw_helper_map( mt->left, depth, mt->left, mt->left->tid);
		//}else
		draw_helper(left(p, base), depth, base, currID);
		for (i = 0; i < (*depth) - 1; i++) printf("  .");

		printf("%3d", _val(p->value)); if (is_marked(p->value)) printf("(x)");

		for (; i < 16; i++) printf("  -");
		printf("<%ld>\n", (long)((long)p - (long)base) / sizeof(struct node));

		//if(p->sright){
		//    struct metadata_struct *mt = ss_lookup(p);
		//    draw_helper_map( mt->right, depth, mt->right, mt->right->tid);
		//}else
		draw_helper(right(p, base), depth, base, currID);
		(*depth)--;
	}
}

void report_all(struct node *p)
{
	int depth = 0;

	draw_helper(p, &depth, p, 0);
}

struct GNode *init_node(int max_node)
{
	struct GNode *first;

#ifdef __PREALLOCGNODES

	/* Atomic element taking from the pool */

	unsigned ctr = atomic_inc(&_poolCtr);

	if (ctr >= MAX_POOLSIZE) {
		printf("Pool exhausted! Exiting...\n");
		exit(27);
	}

	first = &_pool[ctr - 1];

#else

	uintptr_t start = (uintptr_t)calloc(sizeof(struct GNode) + (sizeof(struct node) * (max_node)) +  (sizeof(uintptr_t) * (max_node)), sizeof(char));

	first = (struct GNode*)start;

	start = start + sizeof(struct GNode);

	first->a = (struct node*)start;

	start = start + (sizeof(struct node) * max_node);

	first->b = (uintptr_t*)(start);

#endif

	first->rev = 0;

	gbst_lock_init(&first->lock);

	return first;
}

struct GNode *init_leaf(int max_node)
{
	struct GNode *first = init_node(max_node);

	first->isleaf = 1;
	return first;
}

void out(void *val)
{
	long vals = (long)val;

	printf("%ld ", vals);
}

void printStat(struct global *universe)
{
#ifdef __STATS
	printf("\nNode Count:%d, Node Count(MAX): %d\n", universe->nodecnt, universe->maxcnt);
	printf("\nRebalance (Insert) Done: %d, Rebalance (Delete) Done: %d\n",
	        universe->rebalance_done_ins, universe->rebalance_done_del);
	printf("Insert Count:%d, Delete Count:%d, Failed Insert:%d, Failed Delete:%d \n",
	       universe->count_ins, universe->count_del, universe->failed_ins, universe->failed_del);
#endif
#ifdef __PREALLOCGNODES
	printf("Used nodepool: %d\n", _poolCtr);
#endif

#ifdef __WAIT_COUNT
	printf("Number of waits:%d\n", __wait);
#endif
}


/* integer array printing function */
void print_int_array(const _NODETYPE *array, size_t len)
{
	size_t i;

	printf("\n");
	for (i = 0; i < len; i++)
		printf("%d | ", array[i]);

	printf("\n num: %lu", len);

	printf("\n");
}

/* integer array printing function */
void check_int_array(const int *array, size_t len, struct node *tree)
{
	size_t i;

	int max = 0;

	for (i = 0; i < len; i++) {
		if (array[i] > max) {
			max = array[i];
		} else {
			report_all(tree);
			printf("array error! %d < %d ", array[i], max);
			exit(11);
		}
	}
}


struct node *smart_it_search_lo(struct node *p, void *base, _NODETYPE val)
{
	struct node *last_node = p;

	while (p && p->value != EMPTY) {
		last_node = p;
		if (val < _val(p->value))
			p = left(p, base);
		else
			p = right(p, base);
	}
	return last_node;
}

unsigned sum_node(struct node *p, void *base)
{
	if (!p || p->value == EMPTY)
		return 0;
	else
		return sum_node(left(p, base), base) + sum_node(right(p, base), base) + 1;
}


struct node *smart_it_search(struct node *p, void *base, _NODETYPE val)
{
	struct node *last_node = p;

	while (p && p->value != EMPTY) {
		last_node = p;
		if (val < _val(p->value))
			p = left(p, base);
		else if (_val(p->value) < val)
			p = right(p, base);
		else
			return p;
	}
	return last_node;
}


struct GNode *smart_btree_search_lo(struct GNode *start, _NODETYPE val, int max_depth)
{
	int bits, depth;
	struct node *p;

	//int hop = 0;
	while (!start->isleaf) {
		p = start->a;
		void *base = p;
		uintptr_t *link = start->b;
		bits = 0;
		depth = 0;
		//hop++;
		while (p && p->value != EMPTY) {
			depth++;
			bits <<= 1;
			if (val < _val(p->value)) {
				p = left(p, base);
			} else {
				p = right(p, base);
				bits++;
			}
		}
		bits >>= 1;
		bits <<= ((max_depth - depth));

		//printf("finding %d: bits:%d, depth:%d\n", val, bits, depth);
		if (start->high_key && start->high_key <= val) {
			start = start->sibling;
		} else {
			if (link[bits])
				start = (struct GNode*) link[bits];
			else
				return start;
		}
	}
	//printf( "Hops(val): %d(%d)\n", hop + 1, val);
	return start;
}

int scannode(_NODETYPE key, uintptr_t *temp, int max_depth)
{
	unsigned rev;
	struct GNode *A = (struct GNode*) *temp;
	int bits = 0;
	int depth = 0;

	struct node *p = A->a;
	void *base = p;
	uintptr_t *link = A->b;

	if ((rev = A->rev) & 1) {
		*temp = 0; return 0;
	}

	while (p && p->value != EMPTY) {
		depth++;
		bits <<= 1;
		if (key < _val(p->value)) {
			p = left(p, base);
		} else {
			p = right(p, base);
			bits++;
		}
	}

	bits >>= 1;
	bits <<= ((max_depth - depth));

	/* Follow next_right if high_key is less than searched value*/
	if (A->high_key > 0 && A->high_key <= key) {
		*temp = (uintptr_t) A->sibling;
		if (A->rev - rev) {
			*temp = 0; return 0;
		}
		return 1;
	} else {
		*temp = link[bits];
		if (A->rev - rev) {
			*temp = 0; return 0;
		}
		return 0;
	}
}

struct GNode *move_right_old(_NODETYPE key, struct GNode *t, int max_depth)
{
	struct GNode *current = t;

	if (current->high_key > 0) {
		while (scannode(key, (uintptr_t*)&t, max_depth)) {
			gbst_lock(&t->lock);
			gbst_unlock(&current->lock);
			current = t;
		}
	}
	return current;
}

struct GNode *move_right(_NODETYPE key, struct GNode *t)
{
	struct GNode *current = t;

	while (current->high_key > 0 && current->high_key <= key) {
		t = current->sibling;
		gbst_lock(&t->lock);
		gbst_unlock(&current->lock);
		current = t;
	}

	return current;
}


struct GNode *smart_scannode_lo(struct GNode *start, _NODETYPE val, int max_depth, struct stack *stk)
{
	int bits, depth;
	struct node *p;

	stk->num = 0;

	while (start && !start->isleaf) {
		p = start->a;
		void *base = p;
		uintptr_t *link = start->b;
		bits = 0;
		depth = 0;

		while (p && p->value != EMPTY) {
			depth++;
			bits <<= 1;
			if (val < _val(p->value)) {
				p = left(p, base);
			} else {
				p = right(p, base);
				bits++;
			}
		}
		bits >>= 1;
		bits <<= ((max_depth - depth));

		if (start->sibling && start->high_key <= val) {
			start = start->sibling;
		} else {
			if (link[bits]) {
				push(stk, start);
				start = (struct GNode*) link[bits];
			} else {
				return start;
			}
		}
	}

	return start;
}


struct node *detailed_node_search(struct node *p, void *base, _NODETYPE val, int *bits, int *depth)
{
	struct node *last_node = p;

	*bits = 0;
	*depth = 0;
	while (p && p->value != EMPTY) {
		(*depth)++;
		(*bits) <<= 1;
		last_node = p;
		if (val < _val(p->value)) {
			p = left(p, base);
		} else {
			p = right(p, base);
			(*bits)++;
		}
	}
	(*bits) >>= 1;
	return last_node;
}

void smart_fill_val_parent_lo(struct node *p, void *base, _NODETYPE *buf, uintptr_t *link, uintptr_t *linkbuf, int l, int r, int rmax, unsigned max_depth)
{
	unsigned mid = 0;
	_NODETYPE val;
	int origin_bits = 0, NEW_bits = 0, moved_bits = 0;
	int depth = 0;
	struct GNode *child;

	if (r < l || l == rmax) return;

	mid = ((unsigned int)l + (unsigned int)r) >> 1;

	val = buf[mid];
	child = (struct GNode*) linkbuf[mid];

	struct node *temp = detailed_node_search(p, p, val, &origin_bits, &depth);

	struct node *ln = left(temp, base);
	struct node *rn = right(temp, base);

	if (!ln || !rn) {
		report_all((struct node *) base);
		printf("Failed filling in rebalance (parent): %d (%d)\n", val, rmax);
		print_int_array(buf, rmax);
		exit(0);
	}

	NEW_bits = origin_bits;
	moved_bits = origin_bits;

	origin_bits = origin_bits << (max_depth - depth);

	depth++;
	if (val < temp->value) {
		ln->value = val;
		rn->value = temp->value;
		NEW_bits <<= 1;
		moved_bits = (moved_bits << 1) + 1;
	} else {
		rn->value = val;
		ln->value = temp->value;
		temp->value = val;
		NEW_bits = (NEW_bits << 1) + 1;
		moved_bits <<= 1;
	}

	//Adjust the links
	moved_bits = moved_bits << (max_depth - depth);
	NEW_bits = NEW_bits << (max_depth - depth);

	//Place the links
	link[moved_bits] = link[origin_bits];
	link[NEW_bits] = (uintptr_t) child;

	smart_fill_val_parent_lo(p, base, buf, link, linkbuf, l, mid - 1, rmax, max_depth);
	smart_fill_val_parent_lo(p, base, buf, link, linkbuf, mid + 1, r, rmax, max_depth);
}


void fill_linkbuf(uintptr_t *link, uintptr_t *linkbuf, unsigned maxlink, unsigned *counter)
{
	unsigned i = 0;
	*counter = 0;
	for (i = 0; i < maxlink; i++) {
		if (link[i] != 0) {
			linkbuf[(*counter)] = link[i];
			*counter = *counter + 1;
		}
	}
}


void smart_fill_val_parent(struct node *p, void *base, _NODETYPE *buf, uintptr_t *b, uintptr_t *buflink, int l, int r, int rmax)
{
	unsigned mid = 0;

	if (r < l || l == rmax) return;

	mid = ((unsigned int)l + (unsigned int)r) >> 1;

	smart_fill_val_parent(left(p, base), base, buf, b, buflink, l, mid - 1, rmax);

	p->value = buf[mid];

	smart_fill_val_parent(right(p, base), base, buf, b, buflink, mid + 1, r, rmax);
}

void smart_fill_val(struct node *p, void *base, _NODETYPE *buf, int l, int r, int rmax)
{
	unsigned mid = 0;

	if (r < l || l == rmax) return;

	mid = ((unsigned int)l + (unsigned int)r) >> 1;

	smart_fill_val(left(p, base), base, buf, l, mid - 1, rmax);

	p->value = buf[mid];

	smart_fill_val(right(p, base), base, buf, mid + 1, r, rmax);
}

void smart_fill_val_lo(struct node *p, void *base, _NODETYPE *buf, int l, int r, int rmax)
{
	unsigned mid = 0;
	_NODETYPE val;

	if (r < l || l == rmax) return;

	mid = ((unsigned int)l + (unsigned int)r) >> 1;

	val = buf[mid];

	struct node *temp = smart_it_search_lo(p, base, val);

	struct node *ln = left(temp, base);
	struct node *rn = right(temp, base);

	if (!ln || !rn) {
		report_all((struct node *)	base);
		printf("Failed filling in rebalance (leaf):%d\n", val);
		print_int_array(buf, rmax);
		exit(0);
	}

	if (val < temp->value) {
		ln->value = val;
		rn->value = temp->value;
	} else {
		rn->value = val;
		ln->value = temp->value;
		temp->value = val;
	}

	smart_fill_val_lo(p, base, buf, l, mid - 1, rmax);
	smart_fill_val_lo(p, base, buf, mid + 1, r, rmax);
}

void fill_buf_parent(struct node *p, void *base, _NODETYPE *array, uintptr_t *linkarray, uintptr_t *links, unsigned *counter)
{
	if (p && p->value != EMPTY) {
		fill_buf_parent(left(p, base), base, array, linkarray, links, counter);
		if (!is_marked(p->value))
			array[(*counter)++] = p->value;
		fill_buf_parent(right(p, base), base, array, linkarray, links, counter);
	}
}

void fill_buf(struct node *p, void *base, _NODETYPE *array, unsigned *counter)
{
	if (p && p->value != EMPTY) {
		fill_buf(left(p, base), base, array, counter);
		if (!is_marked(p->value))
			array[(*counter)++] = p->value;
		fill_buf(right(p, base), base, array, counter);
	}
}


void fill_buf_lo(struct node *p, void *base, _NODETYPE *array, unsigned *counter)
{
	struct node *l = left(p, base);
	struct node *r = right(p, base);

	if (l && r && l->value != EMPTY && r->value != EMPTY) {
		fill_buf_lo(l, base, array, counter);
		fill_buf_lo(r, base, array, counter);
	} else {
		if (!is_marked(p->value))
			array[(*counter)++] = p->value;
	}
}


int deleteNode_lo(struct global *universe, _NODETYPE val)
{
	int success;
	success = 1;

	DEBUG_PRINT("Deleting: %d", val);

	if (!*universe->root) {
		success = 0;   //Invalid
	} else {
		struct stack stk;
		struct GNode *start = smart_scannode_lo(*universe->root, val, universe->max_depth, &stk);

		if (!start) return success;

		gbst_lock(&start->lock);
		start = move_right(val, start);

		struct node *temp = smart_it_search(start->a, start->a, val);

		if (val == temp->value) {
			temp->value = (temp->value | 1 << 31);
#ifdef __STATS
			atomic_inc(&start->deleted_node);
#endif
		} else {
			success = 0;
		}

		gbst_unlock(&start->lock);
	}
#ifdef __STATS
	if (success) {
		atomic_dec(&universe->nodecnt);
		atomic_inc(&universe->count_del);
	} else {
		atomic_inc(&universe->failed_del);
	}
#endif

#ifdef DEBUG
	if (success)
		DEBUG_PRINT("--Deleted: %d\n", val);
	else
		DEBUG_PRINT("--Failed Deleting: %d, reason %d\n", val, success);

#endif

	return success;
}


void *getData_lo(struct global *universe, _NODETYPE val)
{
	if (!*universe->root) return 0;

	struct GNode *dt = smart_btree_search_lo(*universe->root, val, universe->max_depth);

	struct node *temp = smart_it_search_lo(dt->a, dt->a, val);

	if (temp)
		if (temp->value == val)
			return (void*) dt->b[0];

	return 0;
}


int searchNode_lo(struct global *universe, _NODETYPE val)
{
	struct GNode *dt = smart_btree_search_lo(*universe->root, val, universe->max_depth);

	struct node *temp = smart_it_search(dt->a, dt->a, val);

	if (temp)
		if (temp->value == val)
			return 1;

	return 0;
}

void smart_fill_val_inc(struct node *p, void *base, _NODETYPE *buf, int l, int r)
{
	unsigned mid = 0;

	if (p == 0) return;

	if (l > r) {
		smart_fill_val_inc(left(p, base), base, buf, l, r);

		p->value = EMPTY;

		smart_fill_val_inc(right(p, base), base, buf, l, r);
	} else {
		mid = ((unsigned int)l + (unsigned int)r) >> 1;

		smart_fill_val_inc(left(p, base), base, buf, l, mid - 1);

		p->value = buf[mid];

		smart_fill_val_inc(right(p, base), base, buf, mid + 1, r);
	}
}


int smart_rec_insert(struct GNode *current, struct node *p, void *base, int *depth, int max_depth, _NODETYPE* success, _NODETYPE val)
{
	if (!p) {
		//printf("no space %d\n", current->count_node);
		return 1;
	}

	if (p->value == EMPTY) {
		p->value = val;
		current->count_node++;
		(*success) = val;

		DEBUG_PRINT( "done! %d\n", current->count_node);
		return 0;
	}

	if (val == _val(p->value)) {
		DEBUG_PRINT( "exists!\n");
		(*success) = val;
		return 0;
	} else {

		int ret = 0;

		if (val < _val(p->value)) {
			ret = smart_rec_insert(current, left(p, base), base, depth, max_depth, success, val);

			if (ret > 0) {
				ret += sum_node(right(p, base), base);
				ret ++;

				(*depth)++;

				//float calc = ((float)(1 << *depth) - 1) * (1.0 - ((((float)*depth) / max_depth) * 0.6));

				int calc = (int)((1<<(*depth)) * 0.5) - 1;

				if (ret - 1 < calc && *depth > 2) {
					DEBUG_PRINT(" Need rebalance! %d %d %d\n", ret, calc, *depth);

					unsigned count = 0;

#ifdef __PREALLOCGNODES
					_NODETYPE buf[__PREALLOCGNODES];
#else
					int max_node = (1 << (*depth )) - 1;
					_NODETYPE buf[max_node];
#endif
					fill_buf(p, base, buf, &count);

					DEBUG_PRINT("Count: %d\n", count);

					atomic_inc(&current->rev);

					smart_fill_val_inc(p, base, buf, 0, count - 1);

					atomic_inc(&current->rev);
					
					// Re-Insert
					ret = smart_rec_insert(current, p, base, depth, max_depth, success, val);

					return 0;
				} else {
					return ret;
				}
			} else {
				return ret;
			}
		} else {
			ret = smart_rec_insert(current, right(p, base), base, depth, max_depth, success, val);

			if (ret > 0) {
				ret += sum_node(left(p, base), base);
				ret ++;

				(*depth)++;

				//float calc = ((float)(1 << *depth) - 1) * (1.0 - ((((float)*depth) / max_depth) * 0.6));

				int calc = (int)((1<<(*depth)) * 0.5) - 1;

				if (ret - 1 < calc && *depth > 2) {
					DEBUG_PRINT(" Need rebalance! %d %d %d\n", ret, calc, *depth);

					unsigned count = 0;

#ifdef __PREALLOCGNODES
					_NODETYPE buf[__PREALLOCGNODES];
#else
					int max_node = (1 << (*depth)) - 1;
					_NODETYPE buf[max_node];
#endif
					fill_buf(p, base, buf, &count);

					DEBUG_PRINT("Count: %d\n", count);

					atomic_inc(&current->rev);

					smart_fill_val_inc(p, base, buf, 0, count - 1);

					atomic_inc(&current->rev);

					// Re-Insert
					ret = smart_rec_insert(current, p, base, depth, max_depth, success, val);

					return 0;
				} else {
					return ret;
				}
			} else {
				return ret;
			}
		}
	}
}

void do_insert_leaf(struct GNode * current, int max_depth, _NODETYPE* success, const _NODETYPE key){

	int tempdepth = 0;

	DEBUG_PRINT("Insert %d\n", key);
	smart_rec_insert(current, current->a, current->a, &tempdepth, max_depth, success, key);

#ifdef DEBUG
	if (!(*success==key)) {
		DEBUG_PRINT("Fail leaf!\n");
		exit(12);
	}
#endif

}

void do_insert_node(struct GNode * current, struct GNode * sibling, int max_depth, int max_node, _NODETYPE keyL, _NODETYPE keyR, _NODETYPE *buf, uintptr_t *buflink)
{
	struct node *temp = 0, *templ = 0, *tempr = 0;
	int origin_bits = 0, depth = 0;
	unsigned count = 0, countlink = 0, mid = 0;

	temp = detailed_node_search(current->a, current->a, keyR, &origin_bits, &depth);

	if (left(temp, current->a) == 0) { //REBALANCE
		DEBUG_PRINT("Parent needs rebalance!\n");
		memset(buf, 0, max_node * sizeof(_NODETYPE));
		memset(buflink, 0, max_node * sizeof(void *));

		fill_buf_lo(current->a, current->a, buf, &count);

		fill_linkbuf(current->b, buflink, max_node, &countlink);

		atomic_inc(&current->rev);

		memset(current->a, 0, max_node * sizeof(struct node));
		memset(current->b, 0, max_node * sizeof(void *));

		mid = ((unsigned int)count) >> 1;

		current->a[0].value = buf[mid];
		current->b[0] = buflink[mid];

		smart_fill_val_parent_lo(current->a, current->a, buf, current->b, buflink, 0, mid - 1, count, max_depth);
		smart_fill_val_parent_lo(current->a, current->a, buf, current->b, buflink, mid + 1, count, count, max_depth);

		atomic_inc(&current->rev);

		current->count_node = count * 2;
	#ifdef __STATS
		atomic_inc(&universe->rebalance_done_ins);
	#endif
		origin_bits = 0;
		depth = 0;
		temp = detailed_node_search(current->a, current->a, keyR, &origin_bits, &depth);

#ifdef DEBUG
		if (left(temp, current->a) == 0){
			DEBUG_PRINT("Fail node\n");
			exit(22);
		}
#endif
	}

	templ = left(temp, current->a);
	tempr = right(temp, current->a);
	depth++;

	if (keyR < temp->value) {
		temp->value = keyR;
		templ->value = keyL;
		tempr->value = keyR;
	} else {
		templ->value = temp->value; //Careful!
		temp->value = keyR;
		tempr->value = keyR;
	}

	current->count_node += 2;
	origin_bits = (origin_bits << 1) + 1;

	//Adjust the links
	origin_bits = origin_bits << (max_depth - depth);

	//Place the links
	current->b[origin_bits] = (uintptr_t) sibling;

}

int insert_par(struct global *universe, const _NODETYPE key, void *data)
{
	struct GNode *current= 0, *sibling = 0;
	struct GNode *oldCurrent = 0;
	
	_NODETYPE keyRNext = 0 , keyR = 0;
	_NODETYPE keyLNext = 0, keyL = 0;
	_NODETYPE success = 0;

	int max_node = universe->max_node;
	int max_depth = universe->max_depth;

#ifdef __PREALLOCGNODES
	_NODETYPE buf [__PREALLOCGNODES];
	uintptr_t buflink [__PREALLOCGNODES];
#else
	_NODETYPE buf [max_node];
	uintptr_t buflink [max_node];
#endif
	struct GNode **root = universe->root;

	/* The tree does not exist yet.
	 * Start new tree.
	 */

	if (*root == NULL) {
		//Try lock the global tree
		gbst_lock(&universe->lock);
		if (*root == NULL) {
			*root = init_leaf(max_node);
			(*root)->a->value = key;
			(*root)->count_node = 1;
			gbst_unlock(&universe->lock);

			//report_all ((*root)->a); printf("\n");
			return 1;
		} else {
			//Wait here first
			gbst_unlock(&universe->lock);
		}
	}

	struct stack Nstack;

	Nstack.num = 0;

	current = *root;

	while (!current->isleaf) {
		oldCurrent = current;
		if (!scannode(key, (uintptr_t*)&current, max_depth)) {
			if (!current)
				current = oldCurrent;
			else
				push(&Nstack, oldCurrent);
		}
	}

	gbst_lock(&current->lock);
	current = move_right(key, current);

	while (1) {
		if (current->isleaf) { //LEAF INSERT
			if(current->count_node >= universe->split_thres) {
				//SPLIT LEAF NODE
				DEBUG_PRINT("Split leaf\n");

				unsigned count = 0, split = 0;

				memset(buf, 0, max_node * sizeof(_NODETYPE));

				fill_buf(current->a, current->a, buf, &count);

				atomic_inc(&current->rev);

				memset(current->a, 0, max_node * sizeof(struct node));

				split = (((unsigned int)count) >> 1);

				smart_fill_val(current->a, current->a, buf, 0, split - 1, split);

				// NEW Node
				sibling = init_leaf(max_node);
				atomic_inc(&sibling->rev);

				// Link the siblings
				sibling->high_key = current->high_key;
				sibling->sibling = current->sibling;

				current->high_key = buf[split];
				current->sibling = sibling;

				atomic_inc(&current->rev);

				smart_fill_val(sibling->a, sibling->a, buf, split, count - 1, count);

				atomic_inc(&sibling->rev);

				//Revise counts
				current->count_node = split;
				sibling->count_node = count - split;

				keyR = buf[split];
				keyL = buf[0];

				if(key >= current->high_key)
					do_insert_leaf(current->sibling, max_depth, &success, key);
				else
					do_insert_leaf(current, max_depth, &success, key);

			}else{

				do_insert_leaf(current, max_depth, &success, key);

				gbst_unlock(&current->lock);

				return success;

			}
		} else { //NODE INSERT
			unsigned count = 0, countlink = 0, mid = 0, split = 0;

			if(current->count_node >= universe->split_thres) {
				DEBUG_PRINT("Split parent\n");

				struct GNode *oldSibling = sibling;

				memset(buf, 0, max_node * sizeof(_NODETYPE));
				memset(buflink, 0, max_node * sizeof(void *));

				fill_buf_lo(current->a, current->a, buf, &count);
				fill_linkbuf(current->b, buflink, max_node, &countlink);

				atomic_inc(&current->rev);

				memset(current->a, 0, max_node * sizeof(struct node));
				memset(current->b, 0, max_node * sizeof(uintptr_t *));

				split = (((unsigned int)count) >> 1);

				mid = (((unsigned int)split) >> 1);

				current->a[0].value = buf[mid];
				current->b[0] = buflink[mid];

				smart_fill_val_parent_lo(current->a, current->a, buf, current->b, buflink, 0, mid - 1, split, max_depth);
				smart_fill_val_parent_lo(current->a, current->a, buf, current->b, buflink, mid + 1, split, split, max_depth);

				// NEW Node
				sibling = init_node(max_node);
				atomic_inc(&sibling->rev);

				// Link the siblings
				sibling->high_key = current->high_key;
				sibling->sibling = current->sibling;

				current->high_key = buf[split];
				current->sibling = sibling;

				atomic_inc(&current->rev);

				mid = (((unsigned int)split) >> 1) + split;

				sibling->a->value = buf[mid];
				sibling->b[0] = buflink[mid];

				smart_fill_val_parent_lo(sibling->a, sibling->a, buf, sibling->b, buflink, split, mid - 1, count, max_depth);
				smart_fill_val_parent_lo(sibling->a, sibling->a, buf, sibling->b, buflink, mid + 1, count, count, max_depth);

				atomic_inc(&sibling->rev);

				//Revise counts
				current->count_node = split * 2;
				sibling->count_node = (count - split) * 2;

				keyRNext = buf[split];
				keyLNext = buf[0];
#ifdef __STATS
				atomic_inc(&universe->rebalance_done_del);
#endif
				if(keyR >= current->high_key)
					do_insert_node(current->sibling, oldSibling,  max_depth,  max_node, keyL, keyR, buf, buflink);
				else
					do_insert_node(current, oldSibling,  max_depth,  max_node, keyL, keyR, buf, buflink);

				keyR = keyRNext;
				keyL = keyLNext;
			} else {

				do_insert_node(current, sibling,  max_depth,  max_node, keyL, keyR, buf, buflink);
				gbst_unlock(&current->lock);
				return success;			
			}
		}

		//Now we have to create a NEW root, restrict only to 1 thread
		if (Nstack.num == 0) {
			struct node *templ = 0, *tempr = 0;
			struct GNode *parent = init_node(max_node);
			parent->a->value = keyR;
			parent->b[0] = (uintptr_t) current;

			//Left&right insert
			templ = left(parent->a, parent->a);
			tempr = right(parent->a, parent->a);

			templ->value = keyL; //Careful!
			tempr->value = keyR;

			parent->b[1 << (max_depth - 2)] = (uintptr_t) sibling;

			parent->count_node = 3;
			*root = parent;

			gbst_unlock(&current->lock);
			return success;
		} else {
	
			oldCurrent = current;
		
			current = pop(&Nstack);

			gbst_lock(&current->lock);

			current = move_right(keyR, current);

			gbst_unlock(&oldCurrent->lock);
		}
	}
	return 0;
}

void initial_add(struct global *universe, int num, int range)
{
	int i = 0;

	_NODETYPE j = 0;

	while (i < num) {
		j = (rand() % range) + 1;
		i += insert_par(universe, j, 0);
	}
#ifdef __STATS
	universe->nodecnt = i;
	universe->maxcnt = i;
#endif
}


greenbst_t *greenbst_alloc(int UB)
{
	/* Allocate the universe */
	greenbst_t *universe = (greenbst_t *) malloc(sizeof(struct global));

	universe->root = NULL;

#ifdef __STATS
	universe->failed_ins = 0;
	universe->failed_del = 0;
	universe->count_ins = 0;
	universe->count_del = 0;
	universe->rebalance_done_ins = 0;
	universe->rebalance_done_del = 0;
	universe->nodecnt = 0;                                  // Initial node count
	universe->maxcnt = 0;                                   // Initial MAXnode count
#endif

	universe->max_depth = ceil(log(UB) / log(2));           // Maximum tree depth based on UB

	universe->max_node = (1 << universe->max_depth) - 1;    // Maximum node for the tree

	universe->split_thres = ((universe->max_node + 1) >> 1) - 1;

	printf("Global start!\n");

	gbst_lock_init(&universe->lock);

	printf("Init map start!\n");

	if (!_map)
		init_tree_map(universe, universe->max_node, universe->max_depth);

	universe->root = (struct GNode **) calloc(1, sizeof(struct GNode *));

	unsigned size = 0;
	size += (sizeof(struct GNode));
	size += ((universe->max_node + 1) * sizeof(struct node));
	size += ((universe->max_node + 1) * sizeof(struct GNode *));

	printf("GNode size is: %u bytes\n", size);

#ifdef __USEMUTEX
	printf("GNode lock is using: Mutex\n\n");
#else
	printf("GNode lock is using: Spinlock\n\n");
#endif

#ifdef __WAIT_COUNT
	__wait = 0;
#endif


	return universe;
}

int greenbst_insert(greenbst_t *map, _NODETYPE key, void *data)
{
	_NODETYPE val = key;

	return insert_par(map, val, data);
}

int greenbst_contains(greenbst_t *map, _NODETYPE key)
{
	_NODETYPE val = key;

	return searchNode_lo(map, val);
}

void *greenbst_get(greenbst_t *map, _NODETYPE key)
{
	_NODETYPE val = key;

	return getData_lo(map, val);
}

int greenbst_delete(greenbst_t *map, _NODETYPE key)
{
	_NODETYPE val = key;

	return deleteNode_lo(map, val);
}


struct map *_map = NULL;

#ifdef __PREALLOCGNODES

struct map mapcontent[__PREALLOCGNODES];

//Pool
struct GNode _pool [MAX_POOLSIZE];

//Counter
unsigned _poolCtr = 0;

#else

struct map *mapcontent;

#endif

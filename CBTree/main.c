/*
 main.c

 CBTree - A concurrent B+Tree based on Lehman & Yao paper
 
 by Ibrahim Umar

 Based on the bpt: B+ Tree implementation
 
 bpt copyright is as below:

 */


/*
 *  bpt.c
 */
#define Version "1.13"
/*
 *
 *  bpt:  B+ Tree Implementation
 *  Copyright (C) 2010  Amittai Aviram  http://www.amittai.com
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 *  Author:  Amittai Aviram
 *    http://www.amittai.com
 *    amittai.aviram@gmail.edu or afa13@columbia.edu
 *    Senior Software Engineer
 *    MathWorks, Inc.
 *    3 Apple Hill Drive
 *    Natick, MA 01760
 *  Original Date:  26 June 2010
 *  Last modified: 15 April 2014
 *
 *  This implementation demonstrates the B+ tree data structure
 *  for educational purposes, includin insertion, deletion, search, and display
 *  of the search path, the leaves, or the whole tree.
 *
 *  Must be compiled with a C99-compliant C compiler such as the latest GCC.
 *
 *  Usage:  bpt [order]
 *  where order is an optional argument
 *  (integer MIN_ORDER <= order <= MAX_ORDER)
 *  defined as the maximal number of pointers in any node.
 *
 */

// Uncomment the line below if you are compiling on Windows.
// #define WINDOWS

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

#include <sys/time.h>
#include "common.h"

#include "bench.h"
#include "locks.h"

#ifdef WINDOWS
#define bool char
#define false 0
#define true 1
#endif

#define atomic_inc(P) __sync_add_and_fetch((P), 1)
#define atomic_dec(P) __sync_add_and_fetch((P), -1)

// TYPES.

/* Type representing the record
 * to which a given key refers.
 * In a real B+ tree system, the
 * record would hold data (in a database)
 * or a file (in an operating system)
 * or some other information.
 * Users can rewrite this part of the code
 * to change the type and content
 * of the value field.
 */
typedef struct record {
	int value;
} record;

// GLOBALS.

/* The order determines the maximum and minimum
 * number of entries (keys and pointers) in any
 * node.  Every node has at most order - 1 keys and
 * at least (roughly speaking) half that number.
 * Every leaf has as many pointers to data as keys,
 * and every internal node has one more pointer
 * to a subtree than the number of keys.
 * This global variable is initialized to the
 * default value.
 */

#ifndef DEFAULT_ORDER
#define DEFAULT_ORDER 256
#endif

int order = DEFAULT_ORDER;

// FUNCTION PROTOTYPES.

// Output and utility.
int cut( int length );

// Insertion.

record * make_record(int value);
node * make_node( void );
node * make_leaf( void );
node * start_new_tree(int key, record * pointer);


/*---------------START PARALEL------------------*/

#define STACK_MAX 100


struct Stack {
    struct node* data[STACK_MAX];
    int         size;
};

typedef struct Stack Stack;


void Stack_Init(Stack *S)
{
    S->size = 0;
}

struct node* Stack_Top(Stack *S)
{
    if (S->size == 0) {
        fprintf(stderr, "Error TOP: stack empty\n");
        exit(-1);
    }
    
    return S->data[S->size-1];
}

void Stack_Push(Stack *S, struct node* d)
{
    if (S->size < STACK_MAX)
        S->data[S->size++] = d;
    else{
        fprintf(stderr, "Error: stack full\n");
        exit(EXIT_FAILURE);
    }
}

struct node* Stack_Pop(Stack *S)
{
    if (S->size == 0){
        fprintf(stderr, "Error POP: stack empty\n");
        exit(EXIT_FAILURE);
    }
    else
        return S->data[--(S->size)];
}

int scannode(int key, struct node** temp, int leaf){
    
    int i = 0;

    unsigned rev;
    
    struct node *A = *temp;
    
    rev = A->rev;
    rev = __sync_fetch_and_add(&rev, 1);
    if (rev & 1) {
		*temp = 0; return 0;
	}

    while (i < A->num_keys) {
        if (key >= A->keys[i]) i++;
        else break;
    }
    
    /* Follow next_right if high_key is less than searched value*/
    if(A->high_key > 0 && A->high_key <= key){
        *temp = A->right_link;

        if (A->rev - rev) {
			*temp = 0; return 0;
		}
        return 1;
    }else{
        if(leaf){
            for (i = 0; i < A->num_keys; i++)
                if (A->keys[i] == key) break;
            if (i == A->num_keys)
                *temp = 0;
            else
                *temp = (node *)A->pointers[i];
            return 0;
        }else{
            *temp =  (node *)A->pointers[i];
        }
        if (A->rev - rev) {
			*temp = 0; return 0;
		}
        return 1;
    }
    
}

int search_par(struct node* root, int key)
{
    struct node *current = root, *oldCurrent;
    
    if(root == NULL) return 0;
    
    while (!current->is_leaf) {
        oldCurrent = current;
        if(!scannode(key, &current, 0))
            if(!current)
			    current = oldCurrent;
    }
    
    while ((scannode(key, &current, 1))) {
    
    }
    
    if(current){
        struct record *rec = (struct record*) current;
        if (rec->value == key)
            return 1;
    }
    
    return 0;
}

struct node* move_right(int key, struct node* t)
{
	struct node *current = t;

	while (current->high_key > 0 && current->high_key <= key) {
		t = current->right_link;
		pthread_spin_lock(&t->lock);
		pthread_spin_unlock(&current->lock);
		current = t;
	}

    return current;
}



int insert_par( node ** root, int key, int value ) {
    
	record * pointer;
	
	node *oldCurrent = NULL,*current = NULL, *new_leaf = NULL, *old_leaf = NULL, *child = NULL;

    int * temp_keys;
    void ** temp_pointers;
    int insertion_index, split, i, j;

    
	/* Case: the tree does not exist yet.
	 * Start a new tree.
	 */
    
	if (*root == NULL){
		//printf("Try\n");
        //Try lock the global tree
        		//printf("Proceed\n");
        		pointer = make_record(value);
        		*root = start_new_tree(key, pointer);
        		return 1;
    }
	
    Stack Nstack;
    
    Stack_Init(&Nstack);
    
    current = *root;
    
    while (!current->is_leaf) {
        oldCurrent = current;
        if(!scannode(key, &current, 0))
	    current = oldCurrent;
        else
            Stack_Push(&Nstack, oldCurrent);

    }
 
    pthread_spin_lock(&current->lock);
    current = move_right(key, current);
    
    //Now check whether the value exists
    for (i = 0; i < current->num_keys; i++)
        if (current->keys[i] == key) break;
    
    if (i != current->num_keys){
        pointer = current->pointers[i];
        if (pointer->value == key){
            pthread_spin_unlock(&current->lock);
            return 0;
        }
    }
    
    pointer = make_record(value);
    
    while(1){
        
        if (current->num_keys < order - 1) {
            //insert_into_leaf(current, key, pointer);
            
            if(current->is_leaf){
                atomic_inc(&current->rev);
                insertion_index = 0;
                while (insertion_index < current->num_keys && current->keys[insertion_index] < key)
                    insertion_index++;
                
                for (i = current->num_keys; i > insertion_index; i--) {
                    current->keys[i] = current->keys[i - 1];
                    current->pointers[i] = current->pointers[i - 1];
                }
                current->keys[insertion_index] = key;
                current->pointers[insertion_index] = pointer;
                current->num_keys++;
                atomic_inc(&current->rev);
            }else{
                insertion_index = 0;
                
                while (insertion_index <= current->num_keys &&
                       current->pointers[insertion_index] != old_leaf)
                    insertion_index++;
                
                atomic_inc(&current->rev);

                for (i = current->num_keys; i > insertion_index; i--) {
                    current->pointers[i + 1] = current->pointers[i];
                    current->keys[i] = current->keys[i - 1];
                }
                current->pointers[insertion_index + 1] = pointer;
                current->keys[insertion_index] = key;
                current->num_keys++;
                atomic_inc(&current->rev);

            }
            pthread_spin_unlock(&current->lock);
            return 1;
        } else {  // split
            
            if(current->is_leaf){
                
                new_leaf = make_leaf();
                
                temp_keys = malloc( order * sizeof(int) );
                if (temp_keys == NULL) {
                    perror("Temporary keys array.");
                    exit(EXIT_FAILURE);
                }
                
                temp_pointers = malloc( order * sizeof(void *) );
                if (temp_pointers == NULL) {
                    perror("Temporary pointers array.");
                    exit(EXIT_FAILURE);
                }
                
                insertion_index = 0;
                while (insertion_index < order - 1 && current->keys[insertion_index] < key)
                    insertion_index++;
                
                for (i = 0, j = 0; i < current->num_keys; i++, j++) {
                    if (j == insertion_index) j++;
                    temp_keys[j] = current->keys[i];
                    temp_pointers[j] = current->pointers[i];
                }
                
                temp_keys[insertion_index] = key;
                temp_pointers[insertion_index] = pointer;
                
                current->num_keys = 0;
                
                split = cut(order - 1);
                
                atomic_inc(&current->rev);
                atomic_inc(&new_leaf->rev);

                for (i = 0; i < split; i++) {
                    current->pointers[i] = temp_pointers[i];
                    current->keys[i] = temp_keys[i];
                    current->num_keys++;
                }
                
                for (i = split, j = 0; i < order; i++, j++) {
                    new_leaf->pointers[j] = temp_pointers[i];
                    new_leaf->keys[j] = temp_keys[i];
                    new_leaf->num_keys++;
                }
                
                free(temp_pointers);
                free(temp_keys);
                
                new_leaf->pointers[order - 1] = current->pointers[order - 1];
                current->pointers[order - 1] = new_leaf;
                
                for (i = current->num_keys; i < order - 1; i++)
                    current->pointers[i] = NULL;
                for (i = new_leaf->num_keys; i < order - 1; i++)
                    new_leaf->pointers[i] = NULL;
                
                new_leaf->parent = current->parent;
                
                /* High Keys */
                new_leaf->high_key = current->high_key;
                current->high_key = (new_leaf->keys[0]);
                new_leaf->right_link = current->right_link;
                current->right_link = new_leaf;
                
		atomic_inc(&new_leaf->rev);
                atomic_inc(&current->rev);
                
                old_leaf = current;
                
                pointer = (struct record*) new_leaf;
                key = new_leaf->keys[0];
            }else{
                /* First create a temporary set of keys and pointers
                 * to hold everything in order, including
                 * the new key and pointer, inserted in their
                 * correct places.
                 * Then create a new node and copy half of the
                 * keys and pointers to the old node and
                 * the other half to the new.
                 */
                
                temp_pointers = malloc( (order + 1) * sizeof(node *) );
                if (temp_pointers == NULL) {
                    perror("Temporary pointers array for splitting nodes.");
                    exit(EXIT_FAILURE);
                }
                temp_keys = malloc( order * sizeof(int) );
                if (temp_keys == NULL) {
                    perror("Temporary keys array for splitting nodes.");
                    exit(EXIT_FAILURE);
                }
                
                insertion_index = 0;
                
                while (insertion_index <= current->num_keys &&
                       current->pointers[insertion_index] != old_leaf)
                    insertion_index++;
                
                for (i = 0, j = 0; i < current->num_keys + 1; i++, j++) {
                    if (j == insertion_index  + 1) j++;
                    temp_pointers[j] = current->pointers[i];
                }
                
                for (i = 0, j = 0; i < current->num_keys; i++, j++) {
                    if (j == insertion_index ) j++;
                    temp_keys[j] = current->keys[i];
                }
                
                temp_pointers[insertion_index  + 1] = pointer;
                temp_keys[insertion_index ] = key;

                /* Create the new node and copy
                 * half the keys and pointers to the
                 * old and half to the new.
                 */
                split = cut(order);
                new_leaf = make_node();
                current->num_keys = 0;

                atomic_inc(&current->rev);
atomic_inc(&new_leaf->rev);


                for (i = 0; i < split - 1; i++) {
                    current->pointers[i] = temp_pointers[i];
                    current->keys[i] = temp_keys[i];
                    current->num_keys++;
                }
                current->pointers[i] = temp_pointers[i];
                key = temp_keys[split - 1];
                for (++i, j = 0; i < order; i++, j++) {
                    new_leaf->pointers[j] = temp_pointers[i];
                    new_leaf->keys[j] = temp_keys[i];
                    new_leaf->num_keys++;
                }
                new_leaf->pointers[j] = temp_pointers[i];
                
                free(temp_pointers);
                free(temp_keys);
                new_leaf->parent = current->parent;
                for (i = 0; i <= new_leaf->num_keys; i++) {
                    child = new_leaf->pointers[i];
                    child->parent = new_leaf;
                }
                
                /* High Keys & Links */
                new_leaf->high_key = current->high_key;
                current->high_key = (new_leaf->keys[0]);
                new_leaf->right_link = current->right_link;
                current->right_link = new_leaf;
                atomic_inc(&new_leaf->rev);

                atomic_inc(&current->rev);

                
                /* Insert a new key into the parent of the two
                 * nodes resulting from the split, with
                 * the old node to the left and the new to the right.
                 */
                
                //return insert_into_parent(root, old_node, k_prime, new_node);
                old_leaf = current;
                
                pointer = (struct record*) new_leaf;
            }
            //Now we have to create a new root, restrict only 1 thread
            if(Nstack.size == 0){
				//printf(".\n");
				struct node *tmp = make_node();
				tmp->keys[0] = key;
				tmp->pointers[0] = current;
				tmp->pointers[1] = new_leaf;
				tmp->num_keys++;
				tmp->parent = NULL;
				current->parent = *root;
				new_leaf->parent = *root;
				*root = tmp;
				pthread_spin_unlock(&current->lock);
				pthread_spin_unlock(&new_leaf->lock);
				pthread_spin_unlock(&tmp->lock);
				return 1;
            }else{
				current = Stack_Pop(&Nstack);
            	pthread_spin_lock(&current->lock);
            	move_right(key, current);
				pthread_spin_unlock(&old_leaf->lock);
				pthread_spin_unlock(&new_leaf->lock);
 			}
        }
    }
    return 1;
    
}

/* Master deletion function.
 */
int delete_par(node * root, int key) {
    
	int j = 0, i = 0, num_pointers;

    struct node* current = NULL, *oldCurrent = NULL;
    struct record* pointer = NULL;
    
    current = root;
    
    if(current == NULL) return 0;

    while (!current->is_leaf) {
        oldCurrent = current;
        if(!scannode(key, &current, 0))
           current = oldCurrent;
    }
    
    pthread_spin_lock(&current->lock);
    current = move_right(key, current);
    
    //Now check whether the value exists
    for (j = 0; j < current->num_keys; j++)
        if (current->keys[j] == key) break;
    
    if (j != current->num_keys){
        pointer = current->pointers[j];
        if (pointer->value == key){

            atomic_inc(&current->rev);
                
            // Remove the key and shift other keys accordingly.
            i = 0;
            while (current->keys[i] != key)
                i++;
            for (++i; i < current->num_keys; i++)
                current->keys[i - 1] = current->keys[i];
            
            // Remove the pointer and shift other pointers accordingly.
            // First determine number of pointers.
            num_pointers = current->is_leaf ? current->num_keys : current->num_keys + 1;
            i = 0;
            while (current->pointers[i] != pointer)
                i++;
            for (++i; i < num_pointers; i++)
                current->pointers[i - 1] = current->pointers[i];
            
            
            // One key fewer.
            current->num_keys--;
            
            // Set the other pointers to NULL for tidiness.
            // A leaf uses the last pointer to point to the next leaf.
            if (current->is_leaf)
                for (i = current->num_keys; i < order - 1; i++)
                    current->pointers[i] = NULL;
            else
                for (i = current->num_keys + 1; i < order; i++)
                    current->pointers[i] = NULL;

            atomic_inc(&current->rev);

            pthread_spin_unlock(&current->lock);
            return 1;
        }
    }
    pthread_spin_unlock(&current->lock);
    
	return 0;
}


/*------------------END PARALEL------------------*/



void destroy_tree_nodes(node * root) {
	int i;
	if (root->is_leaf)
		for (i = 0; i < root->num_keys; i++)
			free(root->pointers[i]);
	else
		for (i = 0; i < root->num_keys + 1; i++)
			destroy_tree_nodes(root->pointers[i]);
	free(root->pointers);
	free(root->keys);
	free(root);
}


node * destroy_tree(node * root) {
	destroy_tree_nodes(root);
	return NULL;
}

/* Finds the appropriate place to
 * split a node that is too big into two.
 */
int cut( int length ) {
	if (length % 2 == 0)
		return length/2;
	else
		return length/2 + 1;
}


/* First insertion:
 * start a new tree.
 */
node * start_new_tree(int key, record * pointer) {
    
	node * root = make_leaf();
	root->keys[0] = key;
	root->pointers[0] = pointer;
	root->pointers[order - 1] = NULL;
	root->parent = NULL;
	root->num_keys++;
    
	return root;
}

/* Creates a new record to hold the value
 * to which a key refers.
 */
record * make_record(int value) {
	record * new_record = (record *)malloc(sizeof(record));
	if (new_record == NULL) {
		perror("Record creation.");
		exit(EXIT_FAILURE);
	}
	else {
		new_record->value = value;
        //new_record->mark = 0;
	}
	return new_record;
}


/* Creates a new general node, which can be adapted
 * to serve as either a leaf or an internal node.
 */
node * make_node( void ) {
	node * new_node;
	new_node = malloc(sizeof(node));
	if (new_node == NULL) {
		perror("Node creation.");
		exit(EXIT_FAILURE);
	}
	new_node->keys = malloc( (order - 1) * sizeof(int) );
	if (new_node->keys == NULL) {
		perror("New node keys array.");
		exit(EXIT_FAILURE);
	}
	new_node->pointers = malloc( order * sizeof(void *) );
	if (new_node->pointers == NULL) {
		perror("New node pointers array.");
		exit(EXIT_FAILURE);
	}
	new_node->is_leaf = false;
	new_node->num_keys = 0;
	new_node->parent = NULL;
	new_node->next = NULL;
    new_node->high_key = 0;
    new_node->right_link = NULL;
    new_node->rev = 0;
    
    pthread_spin_init(&new_node->lock, PTHREAD_PROCESS_SHARED);
    pthread_spin_lock(&new_node->lock);
	return new_node;
}

/* Creates a new leaf by creating a node
 * and then adapting it appropriately.
 */
node * make_leaf( void ) {
	node * leaf = make_node();
	leaf->is_leaf = true;
	return leaf;
}


void initial_add (struct node **root, int num, int range) {
    int i = 0, j = 0;

    while(i < num){
        j = (rand()%range) + 1;
        i += insert_par(root, j, j);
    }
}



/*------------------- END BENCHMARK ---------------------*/

#include<unistd.h>

int main( int argc, char ** argv ) {


struct node *root = NULL;

int myopt = 0;

int s, u, n, i, t, r, v;       //Various parameters

float d;

i = 1023;           //default initial element count
t = 1023;           //default triangle size
r = 5000000;        //default range size
u = 10;             //default update rate
s = 0;              //default seed
n = 1;              //default number of thread
d = (float)1/2;     //default density

v = 0;              //default valgrind mode (reduce stats)

fprintf(stderr,"\nConcurrent BTree\n===============\n\n");
if(argc < 2)
fprintf(stderr,"NOTE: No parameters supplied, will continue with defaults\n");
fprintf(stderr,"Use -h switch for help.\n\n");

while( EOF != myopt ) {
    myopt = getopt(argc,argv,"r:n:i:u:s:hb:");
    switch( myopt ) {
            case 'r': r = atoi( optarg ); break;
            case 'n': n = atoi( optarg ); break;
            case 'i': i = atoi( optarg ); break;
            case 'u': u = atoi( optarg ); break;
            case 's': s = atoi( optarg ); break;
            case 'h': fprintf(stderr,"Accepted parameters\n");
            fprintf(stderr,"-r <NUM>    : Range size\n");
            fprintf(stderr,"-u <0..100> : Update ratio. 0 = Only search; 100 = Only updates\n");
            fprintf(stderr,"-i <NUM>    : Initial tree size (inital pre-filled element count)\n");
            fprintf(stderr,"-n <NUM>    : Number of threads\n");
            fprintf(stderr,"-s <NUM>    : Random seed. 0 = using time as seed\n");
            fprintf(stderr,"-h          : This help\n\n");
            fprintf(stderr,"Benchmark output format: \n\"0: range, insert ratio, delete ratio, #threads, attempted insert, attempted delete, attempted search, effective insert, effective delete, effective search, time (in msec)\"\n\n");
            exit(0);
    }
}
fprintf(stderr,"Parameters:\n");
fprintf(stderr,"- Range size r:\t\t %d\n", r);
fprintf(stderr,"- Update rate u:\t %d%% \n", u);
fprintf(stderr,"- Number of threads n:\t %d\n", n);
fprintf(stderr,"- Initial tree size i:\t %d\n", i);
fprintf(stderr,"- Random seed s:\t %d\n", s);

if (s == 0)
srand((int)time(0));
else
srand(s);

    fprintf(stderr, "Node size: %lu bytes\n", sizeof(node) + ((order - 1) * sizeof(int)) + (order * sizeof(void *)) );

void *pointer = make_record(1);
root = start_new_tree(1, pointer);
pthread_spin_unlock(&root->lock);

#if !defined(__TEST)

    initial_add(&root, 1, r);

    if (i) {
	fprintf(stderr, "Now pre-filling %d random elements...\n", i);
	start_prefill(&root, r, u, 1, i);
    }

    start_benchmark(&root, r, u, n, v);

#else
       testseq(&root, 1);
	testpar(&root, u, n, 1);
    
#endif

return 0;
}

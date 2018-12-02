/*
 bench.c

 This is part of the tree library

 Copyright 2015 Ibrahim Umar (UiT the Arctic University of Norway)

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
 */


#define _GNU_SOURCE
#include <sched.h>

#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <sys/time.h>
#include <math.h>
#include  <pthread.h>

#include "barrier.h"
#include "bench.h"


#if !defined(__APPLE__) && !defined(__MACH__)
#define __THREAD_PINNING 1
#endif

#if defined( __SIM)

#include "m5op.h"
#include "m5_mem.h"
#define __THREAD_PINNING 0
#define MAX_ITER 10000

#else

#define MAX_ITER 5000000

#endif

// MAXITER
int MAXITER;

// Flag for stop condition
volatile char stopFlag;

#ifdef __USEPROF
#include "../common/papicounters.h"
long long **thread_local_values;
#endif


#ifdef __ENERGY
#include "../common/power.h"
#endif

pthread_barrier_t bench_barrier;

/* RANDOM GENERATOR */


// RANGE: (1 - r)
static inline int rand_range_re(unsigned int *seed, long r) {
    return (rand_r(seed) % r) + 1;
}


/* simple function for generating random integer for probability, only works on value of integer 1-100% */

#define MAX_POOL 1000

int p_pool[MAX_POOL];

void prepare_randintp(float ins, float del) {
    
    int i,j=0;
    
    //Put insert
    for(i = 0;i < ins * MAX_POOL/100; i++){
        p_pool[j++]=1;
    }
    //Put delete
    for(i = 0; i< del* MAX_POOL/100; i++){
        p_pool[j++]=2;
    }
    //Put search
    for(i = j; i < MAX_POOL; i++){
        p_pool[j++]=3;
    }
    
    /*
     fprintf(stderr,"\n");
     for (i = 0; i < MAX_POOL; i++)
     fprintf(stderr, "%d, ", p_pool[i]);
     fprintf(stderr,"\n");
     */
}

/* Struct for data input/output per-thread */
struct arg_bench {
    unsigned rank;
#if (__THREAD_PINNING == 1)
    unsigned pinOrder;
#endif
    int runType;
	int size;
    unsigned seed;
    unsigned seed2;
    unsigned update;
    unsigned *pool;
    unsigned long max_iter;
    long counter_ins;
    long counter_del;
    long counter_search;
    long counter_ins_s;
    long counter_del_s;
    long counter_search_s;
    long timer;
    long *inputs;
    int *ops;
    data_t root;
};


void* do_bench (void* arguments)
{
    long counter[3]={0}, success[3]={0};
    int val = 0;// last = 0;
    unsigned *pool;
    int ops, ret = 0;
    int b_size;
    long cont = 0;
    long max_iter = 0;
    int runType = 0;

    struct timeval start, end;
    struct arg_bench *args;

    args = (struct arg_bench*) arguments;
    
    data_t root = args->root;
    max_iter = args->max_iter;
    b_size = args->size;
    pool = args->pool;
    runType = args->runType;
    
#if defined BBST && !defined __PREALLOCGNODES
    init_threads(root->max_node, root->nb_thread);
#endif

#ifndef __APPLE__
#if (__THREAD_PINNING == 1)
    pthread_t current_thread = pthread_self();
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(args->pinOrder, &cpuset);
    fprintf(stdout, "Pinning to core %d... %s\n", args->pinOrder, pthread_setaffinity_np(current_thread, sizeof(cpu_set_t), &cpuset)==0?"Success":"Failed");
#endif
#endif

#ifdef BSTTK
  ssalloc_init();

#if GC == 1
  alloc = (ssmem_allocator_t*) malloc(sizeof(ssmem_allocator_t));
  assert(alloc != NULL);
  ssmem_alloc_init_fs_size(alloc, SSMEM_DEFAULT_MEM_SIZE, SSMEM_GC_FREE_SET_SIZE, args->rank);
#endif
	
#endif

#ifdef RCUT
	urcu_register(args->rank);
#endif

#ifdef BWTREE
	root->AssignGCID(args->rank);
#endif

#ifdef ABTREE
	const int threadID = args->rank;
#endif

#ifdef __USEPROF
    struct localcounters *profcnt;

    if(runType == 1) {
        profcnt = prof_prepare(0);
        thread_local_values[args->rank] = profcnt->values;
    }
#endif

    pthread_barrier_wait(&bench_barrier);

    gettimeofday(&start, NULL);

    if(runType == 1) {
#ifdef __USEPROF
        prof_start(profcnt);
#endif
    }
    /* Check the flag once in a while to see when to quit. */
    while(cont < max_iter && !stopFlag){

        /*---------------------------*/

        //--For a completely random values (original)
        ops = pool[rand_range_re(&args->seed, MAX_POOL) - 1];
        val = rand_range_re(&args->seed2, b_size);
#ifdef LFBST
        switch (ops){
            case 1: ret = BENCH_INSERT(&root[args->rank], val); break;
            case 2: ret = BENCH_DELETE(&root[args->rank], val); break;
            case 3: ret = BENCH_SEARCH(&root[args->rank], val); break;
            default: exit(0); break;
        }
#else
	switch (ops){
            case 1: ret = BENCH_INSERT(root, val); break;
            case 2: ret = BENCH_DELETE(root, val); break;
            case 3: ret = BENCH_SEARCH(root, val); break;
            default: exit(0); break;
        }
#endif
        cont++;
        counter[ops-1]++;
        if(ret)
            success[ops-1]++;
        
    }

    if(runType == 1){
#ifdef __USEPROF
        prof_end(profcnt);
#endif
        // Instruct other threads to stop
        stopFlag = 1;
    }

    gettimeofday(&end, NULL);

    args->counter_ins = counter[0];
    args->counter_del = counter[1];
    args->counter_search = counter[2];
    
    args->counter_ins_s = success[0];
    args->counter_del_s = success[1];
    args->counter_search_s = success[2];
	
#ifdef __USEPROF
    if(runType == 1){
        args->timer = (profcnt->end_time - profcnt->start_time)/1000000;
    }else{
        args->timer = (end.tv_sec * 1000 + end.tv_usec / 1000) - (start.tv_sec * 1000 + start.tv_usec / 1000);
    }
#else
    args->timer = (end.tv_sec * 1000 + end.tv_usec / 1000) - (start.tv_sec * 1000 + start.tv_usec / 1000);
#endif

    pthread_exit((void*) arguments);
}


int benchmark(data_t root, int threads, int key_range, float ins, float del, int initial, int isSilent){
    pthread_t *pid;
    long *inputs;
    int *ops;

    pthread_attr_t attr;

    int i, k;

    struct arg_bench *args, *arg;

    struct arg_bench result;

    // Ensure maximum iteration got filled up depending on the type of operation
    if(initial > 0)
        MAXITER = initial;
    else
        MAXITER = MAX_ITER;

    // Init the stop flag
    stopFlag = 0;

    args = (struct arg_bench*) calloc(threads, sizeof(struct arg_bench));

    inputs = (long*) calloc(key_range, sizeof(long));
    ops = (int*) calloc(key_range, sizeof(int));

    prepare_randintp(ins, del);

#if (__THREAD_PINNING == 1)
    long ncores = sysconf( _SC_NPROCESSORS_ONLN );
    int midcores = (int)ncores/4;

    printf("Available cores: %ld\n", ncores);
#endif

    for(i = 0; i< threads; i++){
        arg = &args[i];
#if (__THREAD_PINNING == 1)
        if(i > midcores -1 && threads < ncores)
            arg->pinOrder = i + midcores;
        else
            arg->pinOrder = i;
#endif
	arg->rank = i;
        // For prefill mode, use zero
        if(initial > 0)
            arg->runType = 0;
        else
            arg->runType = 1;

        arg->size = key_range;

        arg->update = ins + del;
        arg->seed = rand();
        arg->seed2 = rand();
        arg->pool = (unsigned *) calloc(MAX_POOL, sizeof(unsigned));

        for (k = 0; k < MAX_POOL; k++)
            arg->pool[k] = p_pool[k];

        arg->counter_ins = 0;
        arg->counter_del = 0;
        arg->counter_search = 0;

        arg->counter_ins_s = 0;
        arg->counter_del_s = 0;
        arg->counter_search_s = 0;

        arg->timer = 0;

        arg->inputs = inputs;
        arg->ops = ops;

        arg->max_iter = ceil(MAXITER / threads);
        arg->root = root;
    }

    pid = (pthread_t*) calloc(threads, sizeof(pthread_t));

    pthread_attr_init( &attr );

#ifdef PTHREAD_CREATE_UNDETACHED
	pthread_attr_setdetachstate( &attr, PTHREAD_CREATE_UNDETACHED );
	printf("Setting detach\n");
#endif
#ifdef PTHREAD_SCOPE_SYSTEM
	pthread_attr_setscope( &attr, PTHREAD_SCOPE_SYSTEM );
	printf("Setting scope\n");
#endif

#ifdef __USEPROF
    struct localcounters *profcnt;
    if(initial == 0) {
        profcnt = prof_prepare(1);
        thread_local_values = (long long **) calloc(threads, sizeof(long long**));
    }
#endif

    pthread_barrier_init(&bench_barrier,NULL,threads);

    if(initial == 0) {
        fprintf(stderr, "\nStarting benchmark...\n");
    }else{
        fprintf(stderr, "\nInitial insert with %d iterations and %f %% insert proportion and %d threads\n", initial, ins, threads);
    }

    if(initial == 0) {
#ifdef __ENERGY
        ENERGY_START();
#endif

#ifdef __SIM
        map_m5_mem();
        m5_dumpreset_stats(0,0);
#endif
    }

    for (i = 0; i<threads; i++)
        pthread_create (&pid[i], &attr, &do_bench, &args[i]);

    for (i = 0; i<threads; i++)
        pthread_join (pid[i], NULL);

    if(initial == 0) {
#ifdef __SIM
        m5_dumpreset_stats(0,0);
#endif

#ifdef __ENERGY
        ENERGY_END();
#endif
    }

    pthread_attr_destroy( &attr );

    if(!isSilent) {
        fprintf(stderr, "\n0: %d, %0.2f, %0.2f, %d, ", key_range, ins, del, threads);
    }

    result.counter_del = 0;
    result.counter_del_s = 0;
    result.counter_ins = 0;
    result.counter_ins_s = 0;
    result.counter_search = 0;
    result.counter_search_s = 0;
    result.timer = 0;
    
    for(i = 0; i< threads; i++){
        arg = &args[i];
        
        result.counter_del = result.counter_del + arg->counter_del;
        result.counter_del_s = result.counter_del_s + arg->counter_del_s;
        result.counter_ins = result.counter_ins + arg->counter_ins;
        result.counter_ins_s = result.counter_ins_s + arg->counter_ins_s;
        result.counter_search = result.counter_search + arg->counter_search;
        result.counter_search_s = result.counter_search_s + arg->counter_search_s;
        //result.timer = result.timer + arg->timer;
        if(arg->timer > result.timer)
            result.timer = arg->timer;
        
        /*
         fprintf(stderr, "\n(%d): %ld,%ld,%ld,%ld", i, arg->counter_ins, arg->counter_del, arg->counter_search, arg->timer);
         fprintf(stderr, "\n(%d): %ld,%ld,%ld,%ld", i, arg->counter_ins_s, arg->counter_del_s, arg->counter_search_s, arg->timer);
         */

	/* Free all pools */
	free(arg->pool);
    }

    fprintf(stderr, " %ld, %ld, %ld,", result.counter_ins, result.counter_del, result.counter_search);
    fprintf(stderr, " %ld, %ld, %ld, %ld\n", result.counter_ins_s, result.counter_del_s, result.counter_search_s, result.timer);

    if(!isSilent) {
        printf("#OPS,OPERATIONS,%ld\n", result.counter_ins+result.counter_del+result.counter_search);
#ifdef __USEPROF
        printf("#D,TIME,%ld\n", result.timer);
        prof_print_all_threads(threads, thread_local_values);
#endif
    }

    free(pid);
    free(inputs);
    free(ops);
    free(args);

    return 0;
}

void start_prefill(data_t root, int key_range, int updaterate, int num_thread, int initial){

    if(updaterate == 0)
        updaterate = 50;

    float update = (float)updaterate/2;

    // For initial, since number of insert is only percentage, baloon the initial number
    int insertSize = (int) ((100/update) * initial);

#ifdef ABTREE
    const int threadID = 0;
#endif

    // Do single insert first
#ifdef LFBST
    BENCH_INSERT(&root[0], rand()%key_range);
#else
    BENCH_INSERT(root,  rand()%key_range);
#endif

    benchmark(root, num_thread, key_range, update, update, insertSize, 1);

}

void start_benchmark(data_t root, int key_range, int updaterate, int num_thread, int v){
    
    float update = (float)updaterate/2;
    
    benchmark(root, num_thread, key_range, update, update, 0, v);
    
}


struct testpardata{
	int tid;
	int *bulk;
	int nr;
	data_t *untest;
};

void randomize_array(int *values){

	char * pool = (char*) calloc(MAXITER, sizeof(char));

	int cnt = 0;
	int id;

	while (cnt < MAXITER){
		id = rand() % MAXITER;
		if(pool[id] == 0){
			pool[id] = 1;
			values[cnt++] = id + 1;	
		}
	}
	free(pool);
}


void* do_test (void* args){
    
	struct testpardata* arg = (struct testpardata*)args; 
	int myid = arg->tid;
	int nr = arg->nr;
	int* bulk = arg->bulk;
	data_t untest = *arg->untest;
	
    int i;
    
    int start = (MAXITER/nr) * myid;
    int range = (MAXITER/nr);

    int end;

    if( myid == nr-1 ) end = MAXITER;
    else end = start + range;    

    fprintf(stdout, "id:%d, s:%d, r:%d, e:%d\n", myid, start, range, end);

#if defined BBST && !defined __PREALLOCGNODES
    init_threads(untest->max_node, untest->nb_thread);
#endif

#ifdef ABTREE
    const int threadID = arg->tid;
#endif

#ifdef BSTTK
  ssalloc_init();

#if GC == 1
  alloc = (ssmem_allocator_t*) malloc(sizeof(ssmem_allocator_t));
  assert(alloc != NULL);
  ssmem_alloc_init(alloc, SSMEM_DEFAULT_MEM_SIZE, myid);
#endif
	
#if defined(HTICKET)
  init_thread_htlocks(the_cores[d->id]);
#elif defined(CLH)
  init_clh_thread(&clh_local_p);
#endif

#endif
	
#ifdef RCUT
		urcu_register(myid);
#endif

#ifdef BWTREE
	untest->AssignGCID(arg->tid);
#endif

    pthread_barrier_wait(&bench_barrier);
    
    for (i = start; i < end; i++){
#ifdef LFBST
		BENCH_INSERT(&untest[myid], bulk[i]);
#else		
        BENCH_INSERT(untest,  bulk[i]);
#endif
	}
    pthread_exit((void*) args);
}

void testpar(data_t root, int updaterate, int num_thread, int random){
    
    int i, count = 0;
    struct timeval st,ed;
    pthread_barrier_init(&bench_barrier, NULL, num_thread + 1);

    pthread_t pid[num_thread];
    struct testpardata args [num_thread];

    MAXITER = MAX_ITER;

    printf("Concurrent Test:\n\n");
    printf("Inserting %d (%s) elements...\n", MAXITER, (random?"Random":"Increasing"));

	
    int allkey = MAXITER;
    int *bulk = (int*) calloc(MAXITER, sizeof(int));

    for (i = 0; i<num_thread; i++){
		args[i].tid = i;
		args[i].nr = num_thread;
		args[i].bulk = bulk;
		args[i].untest = &root;
    }

    if(random)
	randomize_array(bulk);
    else
    	for(i = 0; i < allkey; i++){
            bulk[i] = 1 + i;
	}
    
    for (i = 0; i<num_thread; i++){
        pthread_create (&pid[i], NULL, &do_test, &args[i]);
    }
    
    pthread_barrier_wait(&bench_barrier);
    
    gettimeofday(&st, NULL);
    
    for (i = 0; i<num_thread; i++)
        pthread_join (pid[i], NULL);
    
    gettimeofday(&ed, NULL);
    
    printf("time : %lu usec\n", (ed.tv_sec - st.tv_sec)*1000000 + ed.tv_usec - st.tv_usec);
    
    //report_all((*untest->root)->a);

#ifdef ABTREE
    const int threadID = 0; 
#endif

    for(i = 0; i < allkey; i++){
#ifdef LFBST
		if(BENCH_SEARCH(&root[0], bulk[i]) < 1){
#else		
		if(BENCH_SEARCH(root, bulk[i]) < 1){
#endif
    		count++;
		}
    }
    fprintf(stderr, "Error searching :%d!\n",count);
    free(bulk);
}

void testseq(data_t root, int random){

  int i, count = 0, seed;
  struct timeval st,ed; 

  int *values;

  MAXITER = MAX_ITER;

  values = (int *) calloc(MAXITER, sizeof(int));
	
  seed = (int) time(NULL);
  srand(seed);

  if(random)
     randomize_array(values);
  else{
     for(i = 0; i < MAXITER; i++)
         values[i] = i+1;
  }
  printf("Sequential Test:\n\n");
  printf("Inserting %d (%s) elements...\n", MAXITER, (random?"Random":"Increasing"));

  gettimeofday(&st, NULL);

#ifdef ABTREE
  const int threadID = 0;
#endif

  for(i = 0; i < MAXITER; i++){
#ifdef LFBST
	BENCH_INSERT(&root[0], values[i]);
#else		
        BENCH_INSERT(root,  values[i]);
#endif
  }

  gettimeofday(&ed, NULL);

  printf("insert time : %lu usec\n", (ed.tv_sec - st.tv_sec)*1000000 + ed.tv_usec - st.tv_usec);
    
  srand(seed);

  gettimeofday(&st, NULL);

  for(i = 0; i < MAXITER; i++){
#ifdef LFBST
		if(BENCH_SEARCH(&root[0], values[i]) < 1){
#else		
		if(BENCH_SEARCH(root, values[i]) < 1){
#endif
            count++;
        }
  }

  gettimeofday(&ed, NULL);    
  printf("search time : %lu usec\n", (ed.tv_sec - st.tv_sec)*1000000 + ed.tv_usec - st.tv_usec);

  fprintf(stderr, "Error searching :%d!\n",count);

#ifdef __TESTSPARSE
  fprintf(stdout, "The size of pool is %d\n", MAXITER * sizeof(int));

  fprintf(stdout, "Press any key!\n");
  getchar();
  /* 10% deletion */
  for(i = 0; i < (MAXITER - (MAXITER * 0.9)); i++){
#ifdef LFBST
	BENCH_DELETE(&root[0], values[i]);
#else
	BENCH_DELETE(root, values[i]);
#endif
      }

  count = 0;

  for(i = 0; i < MAXITER; i++){
#ifdef LFBST
	if(BENCH_SEARCH(&root[0], values[i]) < 1){
#else		
	if(BENCH_SEARCH(root, values[i]) < 1){
#endif
		count++;
	}
  }

  fprintf(stderr, "Values left: %d\n", count);

  free(values);
  fprintf(stdout, "Press any key to continue!\n");
  getchar();
  exit(0);
#endif

  free(values);
}


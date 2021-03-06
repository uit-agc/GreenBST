#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <getopt.h>
#include <time.h>
#include <assert.h>

#include "bwtree.h"
#include "bench.h"


int getVals(TreeType *t1, unsigned x) {
	auto value_set = t1->GetValue(x);
	if(value_set.size()==0) return 0;
	else return 1;
}

int main(int argc, char **argv)
{
	int myopt = 0;

	int s, u, n, i, t, r, v;        //Various parameters

	i = 1023;                       //default initial element count
	t = 4095;                       //default triangle size
	r = 5000000;                    //default range size
	u = 10;                         //default update rate
	s = 0;                          //default seed
	n = 1;                          //default number of thread

	v = 0;                          //default valgrind mode (reduce stats)

	fprintf(stderr, "\nBWTREE\n===============\n\n");
	if (argc < 2)
		fprintf(stderr, "NOTE: No parameters supplied, will continue with defaults\n");
	fprintf(stderr, "Use -h switch for help.\n\n");

	while (EOF != myopt) {
		myopt = getopt(argc, argv, "r:n:i:u:s:hb:");
		switch (myopt) {
		case 'r': r = atoi(optarg); break;
		case 'n': n = atoi(optarg); break;
		case 'i': i = atoi(optarg); break;
		case 'u': u = atoi(optarg); break;
		case 's': s = atoi(optarg); break;
		case 'h': fprintf(stderr, "Accepted parameters\n");
			fprintf(stderr, "-r <NUM>    : Range size\n");
			fprintf(stderr, "-u <0..100> : Update ratio. 0 = Only search; 100 = Only updates\n");
			fprintf(stderr, "-i <NUM>    : Initial tree size (inital pre-filled element count)\n");
			fprintf(stderr, "-n <NUM>    : Number of threads\n");
			fprintf(stderr, "-s <NUM>    : Random seed. 0 = using time as seed\n");
			fprintf(stderr, "-h          : This help\n\n");
			fprintf(stderr, "Benchmark output format: \n\"0: range, insert ratio, delete ratio, #threads, attempted insert, attempted delete, attempted search, effective insert, effective delete, effective search, time (in msec)\"\n\n");
			exit(0);
		}
	}
	fprintf(stderr, "Parameters:\n");
	fprintf(stderr, "- Range size r:\t\t %d\n", r);
	fprintf(stderr, "- Update rate u:\t %d%% \n", u);
	fprintf(stderr, "- Number of threads n:\t %d\n", n);
	fprintf(stderr, "- Initial tree size i:\t %d\n", i);
	fprintf(stderr, "- Random seed s:\t %d\n", s);

	if (s == 0)
		srand((int)time(0));
	else
		srand(s);

	// Init BWTREE
 	data_t t1 = new TreeType{true, KeyComparator{1}, KeyEqualityChecker{1}};

	t1->UpdateThreadLocal(n);
	t1->AssignGCID(0);

	print_flag = false;
#if !defined(__TEST)
	if (i) {
		fprintf(stderr, "Now pre-filling %d random elements...\n", i);
		start_prefill(t1, r, u, n, i);
	}

	fprintf(stderr, "Finished init a BWTREE with initial %d members\n", i);
	fflush(stderr);

	start_benchmark(t1, r, u, n, v);
#else
	testseq(t1, 1);
	testpar(t1, u, n, 1);

#endif
	exit(0);
}

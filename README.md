To run the benchmark:

1. Please make sure you have the Intel PCM library and PAPI library installed.
2. Edit the ```common/common.mk``` to modify the the Intel PCM and PAPI libraries location.
3. Go to ```bench/``` directory.
4. Run ```make-bins.sh```. This should create the trees' binaries.
5. Run ```runtest.{ARCH}.sh``` to run the benchmark.
6. You will find the benchmark results (in CSV format) inside the ```bench/results/``` directory

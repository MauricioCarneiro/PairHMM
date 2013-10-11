all: base omp avx

profilings:	base-profiling omp-profiling avx-profiling

base-profiling:	.FORCE
	g++ pairhmm-1-base.cc -O2 -W -Wall -fopenmp -pedantic -o pairhmm-1-base -lgomp -lm -pg

omp-profiling:	.FORCE
	g++ pairhmm-2-omp.cc -O2 -W -Wall -fopenmp -pedantic -o pairhmm-2-omp -lgomp -lm -pg

avx-profiling:	.FORCE
	g++ pairhmm-3-avx.cc -O2 -W -Wall -fopenmp -pedantic -o pairhmm-3-avx -lgomp -lm -msse4 -pg

base:	.FORCE
	g++ pairhmm-1-base.cc -O2 -W -Wall -fopenmp -pedantic -o pairhmm-1-base -lgomp -lm

omp:	.FORCE
	g++ pairhmm-2-omp.cc -O2 -W -Wall -fopenmp -pedantic -o pairhmm-2-omp -lgomp -lm

avx:	.FORCE
	g++ pairhmm-3-avx.cc -O2 -W -Wall -fopenmp -pedantic -o pairhmm-3-avx -lgomp -lm -msse4

clean:
	rm -f pairhmm-1-base pairhmm-2-omp pairhmm-3-avx *.o

.FORCE:


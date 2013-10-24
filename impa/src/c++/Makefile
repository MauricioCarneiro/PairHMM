all: base omp hybrid dev

profilings:	base-profiling omp-profiling hybrid-profiling dev-profiling

base-profiling:	.FORCE
	g++ pairhmm-1-base.cc -O2 -W -Wall -fopenmp -pedantic -o pairhmm-1-base -lgomp -lm -pg

omp-profiling:	.FORCE
	g++ pairhmm-2-omp.cc -O2 -W -Wall -fopenmp -pedantic -o pairhmm-2-omp -lgomp -lm -pg

hybrid-profiling:	.FORCE
	g++ pairhmm-3-hybrid-float-double.cc -O2 -W -Wall -fopenmp -pedantic -o pairhmm-3-hybrid-float-double -lgomp -lm -pg

dev-profiling:	.FORCE
	g++ pairhmm-dev.cc -O2 -W -Wall -fopenmp -pedantic -o pairhmm-dev -lgomp -lm -pg

base:	.FORCE
	g++ pairhmm-1-base.cc -O2 -W -Wall -fopenmp -pedantic -o pairhmm-1-base -lgomp -lm

omp:	.FORCE
	g++ pairhmm-2-omp.cc -O2 -W -Wall -fopenmp -pedantic -o pairhmm-2-omp -lgomp -lm

hybrid:	.FORCE
	g++ pairhmm-3-hybrid-float-double.cc -O2 -W -Wall -fopenmp -pedantic -o pairhmm-3-hybrid-float-double -lgomp -lm

dev:	.FORCE
	g++ pairhmm-dev.cc -O2 -W -Wall -fopenmp -pedantic -o pairhmm-dev -lgomp -lm

clean:
	rm -f pairhmm-1-base pairhmm-2-omp pairhmm-dev pairhmm-3-hybrid-float-double

.FORCE:


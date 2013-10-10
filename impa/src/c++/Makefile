all: base omp

base:	.FORCE
	g++ pairhmm-base.cc -O2 -W -Wall -fopenmp -pedantic -o pairhmm-base -lgomp -lm

omp:	.FORCE
	g++ pairhmm-omp.cc -O2 -W -Wall -fopenmp -pedantic -o pairhmm-omp -lgomp -lm

clean:
	rm -f pairhmm-base pairhmm-omp *.o

.FORCE:


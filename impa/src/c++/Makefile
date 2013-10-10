all: base

base:	.FORCE
	g++ pairhmm-base.cc -O2 -W -Wall -fopenmp -pedantic -o pairhmm-base -lgomp -lm

clean:
	rm -f pairhmm-base *.o

.FORCE:


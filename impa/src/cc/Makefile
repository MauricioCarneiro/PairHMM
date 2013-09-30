all:	cc

cc:	.FORCE
	g++ pairhmm.cc common.c -O2 -W -Wall -fopenmp -pedantic -o pairhmm-cc -lgomp -lm

clean:
	rm -f pairhmm-cc *.o

.FORCE:


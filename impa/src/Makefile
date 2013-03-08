all:	cuda c

cuda:	.FORCE
	nvcc pairhmm.cu common.cu -O2 -arch=sm_20 -o pairhmm-cuda -lm

c:	.FORCE
	gcc pairhmm.c common.c -O2 -W -Wall -fopenmp -pedantic -o pairhmm-c -lgomp -lm

clean:
	rm -f pairhmm-c pairhmm-cuda *.o

.FORCE:


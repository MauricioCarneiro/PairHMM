LIBS=/usr/local/cuda-5.0/lib64
INCS=/usr/local/cuda-5.0/include

all:	c

c:	.FORCE
	ln -fs ../common.c common.c
	ln -fs ../common.h common.h
	gcc pairhmm-smalltest.c common.c -O2 -W -Wall -fopenmp -pedantic -o pairhmm-smalltest -lgomp -lm
	rm -f common.c common.h

clean:
	rm -f pairhmm-smalltest common.c common.h

.FORCE:


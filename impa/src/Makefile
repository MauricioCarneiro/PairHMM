LIBS=/usr/local/cuda-5.0/lib64
INCS=/usr/local/cuda-5.0/include

all:	cuda c

cuda:	.FORCE
	ln -fs common.c common.cu
	nvcc --ptxas-options=-v -L$(LIBS) -I$(INCS) pairhmm.cu common.cu -O2 -arch=sm_20 -o pairhmm-cuda -lm

c:	.FORCE
	gcc pairhmm.c common.c -O2 -W -Wall -fopenmp -pedantic -o pairhmm-c -lgomp -lm

clean:
	rm -f pairhmm-c pairhmm-cuda *.o common.cu

.FORCE:


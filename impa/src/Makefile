CUDA_VERSION?=5.0
CUDA=/usr/local/cuda-$(CUDA_VERSION)
LIBS=$(CUDA)/lib64
INCS=$(CUDA)/include
NVCC=$(CUDA)/bin/nvcc

all:	cuda cc

cuda:	.FORCE
	ln -fs common.c common.cu
	$(NVCC) --ptxas-options=-v -L$(LIBS) -I$(INCS) pairhmm.cu common.cu -O2 -arch=sm_20 -o pairhmm-cuda -lm

cc:	.FORCE
	g++ pairhmm.cc common.c -O2 -W -Wall -fopenmp -pedantic -o pairhmm-cc -lgomp -lm

clean:
	rm -f pairhmm-cuda pairhmm-cc *.o common.cu

.FORCE:


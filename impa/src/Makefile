CUDA_VERSION?=5.0
CUDA=/usr/local/cuda-$(CUDA_VERSION)
LIBS=$(CUDA)/lib64
INCS=$(CUDA)/include
NVCC=$(CUDA)/bin/nvcc

all:	cuda cc

cuda:	.FORCE
	ln -fs common.c common.cu
	$(NVCC) --ptxas-options=-v -L$(LIBS) -I$(INCS) pairhmm.cu common.cu -O2 \
		-gencode arch=compute_20,code=sm_20 \
		-gencode arch=compute_30,code=sm_30 \
		-gencode arch=compute_35,code=sm_35 \
		-o pairhmm-cuda -lm

inp:
	ln -fs ../../../../test_data ./inp

cc:	.FORCE
	g++ pairhmm.cc common.c -O2 -W -Wall -fopenmp -pedantic -o pairhmm-cc -lgomp -lm

cleanlocal:
	rm -f pairhmm-cuda pairhmm-cc *.o common.cu inp

clean:
	rm -f pairhmm-cuda pairhmm-cc *.o common.cu inp
	ssh eric@dijkstra 'rm -rf /home/eric/PairHMM/rep/PairHMM/impa/src'
	ssh eric@germain 'rm -rf /home/eric/PairHMM/rep/PairHMM/impa/src'
	ssh eric@bresenham 'rm -rf /home/eric/PairHMM/rep/PairHMM/impa/src'

.FORCE:


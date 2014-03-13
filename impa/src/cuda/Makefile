CUDA_VERSION?=5.5
CUDA=/usr/local/cuda-$(CUDA_VERSION)
LIBS=$(CUDA)/lib64
INCS=$(CUDA)/include
NVCC=$(CUDA)/bin/nvcc

FLAGS=--ptxas-options=-v -L$(LIBS) -I$(INCS) -O2 -gencode arch=compute_20,code=sm_20 -gencode arch=compute_30,code=sm_30 -gencode arch=compute_35,code=sm_35 -lm

#BIN:=pairhmm-1-base pairhmm-dev
BIN:=pairhmm-dev

all:	$(BIN)

pairhmm-1-base:	pairhmm-1-base.cu 
	$(NVCC) -o $@-32 $^ -DNTHREADS=32 $(FLAGS)
	$(NVCC) -o $@-64 $^ -DNTHREADS=64 $(FLAGS)
	$(NVCC) -o $@-96 $^ -DNTHREADS=96 $(FLAGS)

pairhmm-dev:	pairhmm-dev.cu 
	$(NVCC) -o $@-32 $^ -DNTHREADS=32 $(FLAGS)
	$(NVCC) -o $@-64 $^ -DNTHREADS=64 $(FLAGS)
	$(NVCC) -o $@-96 $^ -DNTHREADS=96 $(FLAGS)

clean:
	rm -f $(foreach bin,$(BIN),$(bin)-32)
	rm -f $(foreach bin,$(BIN),$(bin)-64)
	rm -f $(foreach bin,$(BIN),$(bin)-96)
	rm -f *.precision
	rm -f *.err

CUDA_VERSION?=5.5
CUDA=/usr/local/cuda-$(CUDA_VERSION)
LIBS=$(CUDA)/lib64
INCS=$(CUDA)/include
NVCC=$(CUDA)/bin/nvcc

all:	pairhmm32x1 pairhmm32x3 pairhmm64x1 pairhmm64x3 pairhmm96x1 pairhmm96x3

pairhmm32x1:	pairhmm.cu common.cu
	$(NVCC) --ptxas-options=-v -L$(LIBS) -I$(INCS) pairhmm.cu common.cu -O2 \
		-gencode arch=compute_20,code=sm_20 \
		-gencode arch=compute_30,code=sm_30 \
		-gencode arch=compute_35,code=sm_35 \
		-DBLOCKWIDTH=32 -DBLOCKHEIGHT=1 \
		-o pairhmm32x1 -lm

pairhmm32x3:	pairhmm.cu common.cu
	$(NVCC) --ptxas-options=-v -L$(LIBS) -I$(INCS) pairhmm.cu common.cu -O2 \
		-gencode arch=compute_20,code=sm_20 \
		-gencode arch=compute_30,code=sm_30 \
		-gencode arch=compute_35,code=sm_35 \
		-DBLOCKWIDTH=32 -DBLOCKHEIGHT=3 \
		-o pairhmm32x3 -lm

pairhmm64x1:	pairhmm.cu common.cu
	$(NVCC) --ptxas-options=-v -L$(LIBS) -I$(INCS) pairhmm.cu common.cu -O2 \
		-gencode arch=compute_20,code=sm_20 \
		-gencode arch=compute_30,code=sm_30 \
		-gencode arch=compute_35,code=sm_35 \
		-DBLOCKWIDTH=64 -DBLOCKHEIGHT=1 \
		-o pairhmm64x1 -lm

pairhmm64x3:	pairhmm.cu common.cu
	$(NVCC) --ptxas-options=-v -L$(LIBS) -I$(INCS) pairhmm.cu common.cu -O2 \
		-gencode arch=compute_20,code=sm_20 \
		-gencode arch=compute_30,code=sm_30 \
		-gencode arch=compute_35,code=sm_35 \
		-DBLOCKWIDTH=64 -DBLOCKHEIGHT=3 \
		-o pairhmm64x3 -lm

pairhmm96x1:	pairhmm.cu common.cu
	$(NVCC) --ptxas-options=-v -L$(LIBS) -I$(INCS) pairhmm.cu common.cu -O2 \
		-gencode arch=compute_20,code=sm_20 \
		-gencode arch=compute_30,code=sm_30 \
		-gencode arch=compute_35,code=sm_35 \
		-DBLOCKWIDTH=96 -DBLOCKHEIGHT=1 \
		-o pairhmm96x1 -lm

pairhmm96x3:	pairhmm.cu common.cu
	$(NVCC) --ptxas-options=-v -L$(LIBS) -I$(INCS) pairhmm.cu common.cu -O2 \
		-gencode arch=compute_20,code=sm_20 \
		-gencode arch=compute_30,code=sm_30 \
		-gencode arch=compute_35,code=sm_35 \
		-DBLOCKWIDTH=96 -DBLOCKHEIGHT=3 \
		-o pairhmm96x3 -lm

clean:
	rm -f pairhmm32x1 pairhmm32x3 pairhmm64x1 pairhmm64x3 pairhmm96x1 pairhmm96x3

.FORCE:


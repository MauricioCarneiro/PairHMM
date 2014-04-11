#OMPCFLAGS=-fopenmp
#OMPLDFLAGS=-lgomp
USERCFLAGS=

#CFLAGS=-O2 -std=c++11 -W -Wall -march=corei7-avx -Wa,-q            -pedantic $(OMPCFLAGS) -Wno-unknown-pragmas
CFLAGS= -W -Wall -march=corei7 -mfpmath=sse -msse4.2 -pedantic $(OMPCFLAGS) -Wno-unknown-pragmas $(USERFLAGS)

CXXFLAGS=$(CFLAGS)
CC=gcc
CXX=g++
NVCC=nvcc 
CUDAFLAGS= -Xptxas -v,-abi=no -arch=sm_35 $(USERFLAGS)

ifeq ($(DEBUG),1)
   CFLAGS += -g
   CUDAFLAGS += -G -lineinfo
endif

LDFLAGS=-lm $(OMPLDFLAGS) $(CFLAGS)

#BIN:=pairhmm-1-base pairhmm-2-omp pairhmm-3-hybrid-float-double pairhmm-4-hybrid-diagonal pairhmm-5-hybrid-diagonal-homogeneus pairhmm-6-onlythreediags pairhmm-7-presse pairhmm-8-sse #pairhmm-dev

BIN:= pairhmm-cuda

all: $(BIN)

pairhmm-1-base:	pairhmm-1-base.cc input.o
	$(CXX) -o $@ $^ $(LDFLAGS)

pairhmm-2-omp:	pairhmm-2-omp.cc input.o
	$(CXX) -o $@ $^ $(LDFLAGS)

pairhmm-3-hybrid-float-double:	pairhmm-3-hybrid-float-double.cc input.o
	$(CXX) -o $@ $^ $(LDFLAGS)

pairhmm-4-hybrid-diagonal:	pairhmm-4-hybrid-diagonal.cc input.o
	$(CXX) -o $@ $^ $(LDFLAGS)

pairhmm-5-hybrid-diagonal-homogeneus:	pairhmm-5-hybrid-diagonal-homogeneus.cc input.o
	$(CXX) -o $@ $^ $(LDFLAGS)

pairhmm-6-onlythreediags:	pairhmm-6-onlythreediags.cc input.o
	$(CXX) -o $@ $^ $(LDFLAGS)

pairhmm-7-presse:	pairhmm-7-presse.cc input.o
	$(CXX) -o $@ $^ $(LDFLAGS)

pairhmm-8-sse:	pairhmm-8-sse.cc input.o
	$(CXX) -o $@ $^ $(LDFLAGS)

pairhmm-dev: pairhmm-dev.o input.o
	$(CXX) -o $@ $^ $(LDFLAGS)

pairhmm-cuda: pairhmm-cuda.cu input.o
	echo "Building pairhmm-cuda"
	$(NVCC) -g -o $@ $^ $(CUDAFLAGS)

pairhmm-dbg: pairhmm-dbg.o input.o
	$(CXX) -o $@ $^ $(LDFLAGS) -D DEBUG_MODE

bench: pairhmm-dev
	./benchmark.sh

clean:
	rm -f $(BIN) *.o

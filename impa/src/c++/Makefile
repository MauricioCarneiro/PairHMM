#OMPCFLAGS=-fopenmp
#OMPLDFLAGS=-lgomp

#CFLAGS=-O2 -std=c++11 -W -Wall -march=corei7-avx -Wa,-q -pedantic $(OMPCFLAGS) -Wno-unknown-pragmas
CFLAGS=-O2 -std=c++11 -W -Wall -march=corei7 -mfpmath=sse -msse4.2 -pedantic $(OMPCFLAGS) -Wno-unknown-pragmas

CXXFLAGS=$(CFLAGS)
CC=gcc
CXX=g++

LDFLAGS=-lm $(OMPLDFLAGS)

BIN:=pairhmm-6-onlythreediags pairhmm-dev

all: $(BIN)

pairhmm-6-onlythreediags: pairhmm-6-onlythreediags.o
	$(CXX) -o $@ $^ $(LDFLAGS)

pairhmm-dev: pairhmm-dev.o
	$(CXX) -o $@ $^ $(LDFLAGS)

bench: pairhmm-dev
	./benchmark.sh

clean:
	rm -f $(BIN) *.o

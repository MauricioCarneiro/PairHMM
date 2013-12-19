#OMPCFLAGS=-fopenmp
#OMPLDFLAGS=-lgomp

#CFLAGS=-O2 -g -std=c++11 -W -Wall -march=corei7-avx -Wa,-q -pedantic $(OMPCFLAGS) -Wno-unknown-pragmas -Wno-vla
CFLAGS=-O2 -g -std=c++11 -W -Wall -march=corei7 -mfpmath=sse -msse4.2 -pedantic $(OMPCFLAGS) -Wno-unknown-pragmas -Wno-vla

CXXFLAGS=$(CFLAGS)
CC=gcc
CXX=g++
TOGRIND=pairhmm-dev
GRINDIN=performances/benchmark.in
BENCH=pairhmm-dev

LDFLAGS=-lm $(OMPLDFLAGS)

BIN:=pairhmm-dev

all: $(BIN)

pairhmm-6-onlythreediags: pairhmm-6-onlythreediags.o
	$(CXX) -o $@ $^ $(LDFLAGS)

pairhmm-dev: pairhmm-dev.o
	$(CXX) -o $@ $^ $(LDFLAGS)

pairhmm-dev-avx: pairhmm-dev-avx.o
	$(CXX) -o $@ $^ $(LDFLAGS)

pairhmm-7-presse: pairhmm-7-presse.o
	$(CXX) -o $@ $^ $(LDFLAGS)

pairhmm-8-sse: pairhmm-8-sse.o
	$(CXX) -o $@ $^ $(LDFLAGS)

bench: $(BENCH)
	./benchmark.sh $(BENCH)

valgrind: $(TOGRIND)
	valgrind --dsymutil=yes --tool=memcheck -v $(TOGRIND) < $(GRINDIN)

clean:
	rm -f $(BIN) *.o

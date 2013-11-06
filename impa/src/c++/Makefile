FLAGS=-O2 -W -Wall -Wno-vla -fopenmp -pedantic -lgomp -lm
FLAGS_PROFILING=$(FLAGS) -pg
PROFILING=

all: *.cc

*.cc:	.FORCE
ifeq ($(PROFILING), true)
	g++ $@ -o $(basename $@)_profiling.o $(FLAGS_PROFILING)
else
	g++ $@ -o $(basename $@).o $(FLAGS)
endif

clean:	.FORCE
	rm -f *.o

.FORCE:


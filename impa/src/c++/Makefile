FLAGS=-O2 -std=c++11 -W -Wall -Wno-vla -mfpmath=sse -msse4.2 -fopenmp -pedantic -lgomp -lm
#FLAGS=-W -ggdb -Wall -Wno-vla -pedantic -lm
FLAGS_PROFILING=$(FLAGS) -pg
PROFILING=
DEBUG_MODE=

all: *.cc

*.cc:	.FORCE
ifeq ($(PROFILING), true)
	g++ $@ -o $(basename $@)_profiling $(FLAGS_PROFILING)
else
ifeq ($(DEBUG_MODE), true)
	g++ $@ -o $(basename $@) $(FLAGS) -DDEBUG_MODE=1
else
	g++ $@ -o $(basename $@) $(FLAGS)
endif
endif

clean:	.FORCE
	rm pairhmm-1-base pairhmm-2-omp pairhmm-3-hybrid-float-double pairhmm-4-hybrid-diagonal pairhmm-5-hybrid-diagonal-homogeneus pairhmm-6-onlythreediags pairhmm-7-presse pairhmm-dev

.FORCE:


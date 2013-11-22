FLAGS=-O2 -W -Wall -Wno-vla -fopenmp -pedantic -lgomp -lm
FLAGS_PROFILING=$(FLAGS) -pg
PROFILING=

all: *.cc

*.cc:	.FORCE
ifeq ($(PROFILING), true)
	g++ $@ -o $(basename $@)_profiling $(FLAGS_PROFILING)
else
	g++ $@ -o $(basename $@) $(FLAGS)
endif

clean:	.FORCE
	find . -executable -exec rm '{}' \;

.FORCE:


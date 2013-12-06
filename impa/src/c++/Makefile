FLAGS=-O2 -W -Wall -Wno-vla -fopenmp -pedantic -lgomp -lm
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
	find . -executable -exec rm '{}' \;

.FORCE:


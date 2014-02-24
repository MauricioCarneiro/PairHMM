CXXFLAGS=-c -Wall -std=c++11 -g -O0
LDFLAGS=
SOURCES= testcase.cpp testcase_iterator.cpp haplotype.cpp main.cpp diagonals.cpp pairhmm.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=../pairhmm

all: $(SOURCES) $(EXECUTABLE) 

$(EXECUTABLE): $(OBJECTS) 
	$(CXX) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -f *.o ../pairhmm

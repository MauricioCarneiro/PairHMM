CXXFLAGS=-c -Wall -std=c++11
LDFLAGS=
SOURCES= testcase.cpp testcase_iterator.cpp read.cpp haplotype.cpp main.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=../pairhmm

all: $(SOURCES) $(EXECUTABLE) 

$(EXECUTABLE): $(OBJECTS) 
	$(CXX) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -f *.o ../pairhmm

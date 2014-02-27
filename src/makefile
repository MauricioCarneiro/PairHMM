CXXFLAGS=-c -Wall -std=c++11 -g -O0
LDFLAGS=
SOURCES= testcase.cpp testcase_iterator.cpp haplotype.cpp main.cpp diagonals.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=../pairhmm

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(LDFLAGS) $(OBJECTS) -o $@

diagonals.o: diagonals.cpp diagonals.h
haplotype.o: haplotype.cpp haplotype.h utils.h
main.o: main.cpp input_reader.h testcase_iterator.h testcase.h read.h \
 utils.h haplotype.h pairhmm.h pairhmm_impl.h constants.h diagonals.h 
testcase.o: testcase.cpp testcase.h read.h utils.h haplotype.h
testcase_iterator.o: testcase_iterator.cpp testcase_iterator.h testcase.h \
 read.h utils.h haplotype.h

.cpp.o:
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -f *.o ../pairhmm

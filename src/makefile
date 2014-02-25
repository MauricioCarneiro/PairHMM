CXXFLAGS=-c -Wall -std=c++11 -g -O0
LDFLAGS=
SOURCES= testcase.cpp testcase_iterator.cpp haplotype.cpp main.cpp diagonals.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=../pairhmm

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(LDFLAGS) $(OBJECTS) -o $@

test: ../pairhmm
	@../pairhmm < ../test_data/tiny.in | paste - ../test_data/tiny.out  | awk 'BEGIN {m = 0; n = NR} {a = $$2-$$1; a = a < 0? -a: a; if (a > m) {m = a; n = NR} } END { printf("max error %g at line %d\n", m, n)} '

diagonals.o: diagonals.cpp diagonals.h
haplotype.o: haplotype.cpp haplotype.h utils.h
main.o: main.cpp input_reader.h testcase_iterator.h testcase.h read.h \
 utils.h haplotype.h pairhmm.h pairhmm_impl.h constants.h diagonals.h \
 pairhmm.hpp
testcase.o: testcase.cpp testcase.h read.h utils.h haplotype.h
testcase_iterator.o: testcase_iterator.cpp testcase_iterator.h testcase.h \
 read.h utils.h haplotype.h

.cpp.o:
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -f *.o ../pairhmm

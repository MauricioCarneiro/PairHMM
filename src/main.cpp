#include <iostream>
#include "input_reader.h"
#include "testcase_iterator.h"
#include "testcase.h"

using namespace std;

int main (const int argc, char const * const argv[]) {
  InputReader<TestcaseIterator> reader {};
  if (argc == 2) 
    reader.from_file(argv[0]);
  for (auto& testcase : reader) {
    cout << testcase << endl;
  }
  return 0;
}

#include <iostream>
#include <vector>
#include <cmath>

#include "input_reader.h"
#include "testcase_iterator.h"
#include "testcase.h"
#include "pairhmm.h"

using namespace std;


int main (const int argc, char const * const argv[]) {
  auto pairhmm = Pairhmm<float>{};
  InputReader<TestcaseIterator> reader {};
  if (argc == 2) 
    reader.from_file(argv[1]);
  for (auto& testcase : reader) {
    auto results = pairhmm.calculate(testcase);
    for (auto x : results) 
      cout << x << endl;
  }
  return 0;
}

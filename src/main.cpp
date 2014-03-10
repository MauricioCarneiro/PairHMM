#include <iostream>
#include <vector>
#include <cmath>

#include "input_reader.h"
#include "testcase_iterator.h"
#include "testcase.h"
#include "pairhmm.h"
#include "pairhmm_impl.h"
#include "pairhmm_vecimpl.h"
#include "pairhmm_sseimpl.h"
#include "pairhmm_avximpl.h"
#include "chronos.h"

using namespace std;

int main (const int argc, char const * const argv[]) {
  auto pairhmm = Pairhmm<
    //PairhmmImpl<double, Diagonals<double>, Constants<double>>,
    //PairhmmVecImpl<float, Diagonals<float>, Constants<float>, 4>,
    //PairhmmSSEFloatImpl,
    PairhmmAVXFloatImpl,
    //PairhmmImpl<float, Diagonals<float>, Constants<float>>,
    PairhmmImpl<double, Diagonals<double>, Constants<double>>
  >{};
  InputReader<TestcaseIterator> reader {};
  if (argc == 2)
    reader.from_file(argv[1]);
  double computation_time = 0.f;
  Chronos time;
  for (auto& testcase : reader) {
    time.reset();
    auto results = pairhmm.calculate(testcase);
    computation_time += time.elapsed();
    for (auto x : results)
      cout << x << endl;
  }
  std::cerr << "done in " << computation_time << "ms\n";
  return 0;
}

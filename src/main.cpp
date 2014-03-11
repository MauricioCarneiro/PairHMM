#include <iostream>
#include <vector>
#include <cmath>

#ifndef DENORMALS
#include <xmmintrin.h>
#endif

#include "input_reader.h"
#include "testcase_iterator.h"
#include "testcase.h"
#include "pairhmm.h"
#include "pairhmm_scalarimpl.h"
#include "pairhmm_vecimpl.h"
#include "pairhmm_sseimpl.h"
#include "pairhmm_avximpl.h"
#include "chronos.h"

using namespace std;

int main (const int argc, char const * const argv[]) {
#ifndef DENORMALS
    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
    _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
#endif
  auto pairhmm = Pairhmm<    // this should be detected in runtime (we generate all versions and run the best supported one)
    //PairhmmScalarImpl<float>,
    //PairhmmScalarImpl<double>
    //PairhmmSSEFloatImpl,
    //PairhmmScalarImpl<double>
    //PairhmmAVXFloatImpl,
    //PairhmmAVXDoubleImpl
    PairhmmAVXFloat2DiagsImpl,
    PairhmmAVXDouble2DiagsImpl
  >{};
  InputReader<TestcaseIterator> reader {};
  if (argc == 2)
    reader.from_file(argv[1]);
  auto computation_time = 0.f;
  auto timer = Chronos{};
  for (auto& testcase : reader) {
    timer.reset();
    auto results = pairhmm.calculate(testcase);
    computation_time += timer.elapsed();
    for (auto x : results)
      cout << x << endl;
  }
  std::clog << "computation time: " << computation_time << "ms\n";
  return 0;
}

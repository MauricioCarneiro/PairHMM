#include <iostream>
#include <vector>
#include <cmath>
#include <future>

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

constexpr size_t TESTCASES_TO_READ = 100000;

vector<Testcase> read_testcases(InputReader<TestcaseIterator>& reader, const size_t n) {
  auto result = vector<Testcase>{};
  if (n > 0) {
    result.reserve(n);
    auto i = 0u;
    for (const auto& testcase : reader) {
      result.push_back(testcase);
      if (++i == n)
        break;
    }
  }
  return result;
}

template<class F, class A>
vector<double> calculate(const Pairhmm<F,A>& phmm, Testcase testcase) {
  return phmm.calculate(testcase);
}

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
  auto computation_time = 0.f;
  auto timer = Chronos{};
  InputReader<TestcaseIterator> reader {};
  if (argc == 2)
    reader.from_file(argv[1]);

  while (!reader.eof()) {
    clog << "reading " << TESTCASES_TO_READ << " testcases..." << endl;
    const auto testcases = read_testcases(reader, TESTCASES_TO_READ);
    auto results = vector<future<vector<double>>>{}; // inefficient, rewrite this as a class
    clog << "computing... " << endl;
    timer.reset();
    for (const auto& testcase : testcases) 
      results.push_back(async(launch::async, &Pairhmm<PairhmmAVXFloat2DiagsImpl, PairhmmAVXDouble2DiagsImpl>::calculate, pairhmm, testcase));
    computation_time += timer.elapsed();
    clog << "results generated... " << endl;
    for (auto i = 0u; i != results.size(); ++i)
      for(const auto& x : results[i].get()) 
        cout << x << endl;
  }
  std::clog << "computation time: " << computation_time << "ms\n";
  return 0;
}

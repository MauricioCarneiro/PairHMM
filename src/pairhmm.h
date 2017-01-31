#ifndef PAIRHMM_H
#define PAIRHMM_H

#include <vector>

#include "pairhmm_impl.h"
#include "testcase.h"

template <typename FAST, typename CORRECT>
class Pairhmm {
 private:
  FAST fast;
  CORRECT correct;
  double resow_time = 0.f;
 public:
  std::vector<double> calculate(const Testcase& testcase) {
    auto results = fast.calculate(testcase); // calculate all the testcases using the float precision pairhmm implementation
    correct.recalculate(testcase, &results[0]);  // recalculate only the ones that failed with double precision (results = DOUBLE_MIN for failed ones)
    return results;
  }
  std::vector<double> reap() {
    Chronos time; time.reset();
    auto results = fast.reap();
    time.reset();
    correct.rereap(results);
    std::cerr << "DP computation: " << time.elapsed() << " ms." << std::endl;
    std::cerr << "resow time: " << resow_time << " ms." << std::endl;
    //correct.recalculate(testcase, results);
    return results;
  }
  void sow(const Testcase& testcase) {
    Chronos time;
    fast.sow(testcase);
    time.reset();
    correct.resow(testcase);
    resow_time += time.elapsed();
  }
};

#endif

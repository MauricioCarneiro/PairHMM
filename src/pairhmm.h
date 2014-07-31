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
 public:
  std::vector<double> calculate(const Testcase& testcase) {
    auto results = fast.calculate(testcase); // calculate all the testcases using the float precision pairhmm implementation
    correct.recalculate(testcase, &results[0]);  // recalculate only the ones that failed with double precision (results = DOUBLE_MIN for failed ones)
    return results;
  }
  std::vector<double> reap() {
    auto results = fast.reap();
    //TODO 
    //     - move testcase data from fast to correct
    //     - pad read and haplotype for cases which underflowed
    //     - call correct.compute_full_prob for each
    correct.rereap(results);
    //correct.recalculate(testcase, results);
    return results;
  }
  void sow(const Testcase& testcase) {
    fast.sow(testcase);
    correct.resow(testcase);
    //TODO maybe call correct.sow for each testcase so the data will
    //     be in correct when it's time to recalc and reap
  }
};

#endif

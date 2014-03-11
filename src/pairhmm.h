#ifndef PAIRHMM_H
#define PAIRHMM_H

#include <vector>

#include "pairhmm.h"
#include "testcase.h"

template <typename FAST, typename CORRECT>
class Pairhmm {
 public:
  std::vector<double> calculate(const Testcase& testcase) {
    auto results = fast.calculate(testcase);  // calculate all the testcases using the float precision pairhmm implementation
    correct.recalculate(testcase, results, fast.max_original_read_length(), fast.max_padded_read_length()); // recalculate only the ones that failed with double precision (results = DOUBLE_MIN for failed ones)
    return results;
  }

 private:
  FAST fast;
  CORRECT correct;
};

#endif

#ifndef __PAIRHMM__
#define __PAIRHMM__

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
    auto results = fast.calculate(testcase);  // calculate all the testcases using the float precision pairhmm implementation
    correct.recalculate(testcase, results, fast.max_original_read_length, fast.max_padded_read_length, fast.get_padded_haplotypes()); // recalculate only the ones that failed with double precision (results = DOUBLE_MIN for failed ones)
    return results;
  }
};

#endif

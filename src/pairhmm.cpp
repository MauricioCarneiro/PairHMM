#include <vector>
#include <limits>

#include "pairhmm.h"
#include "testcase.h"

using namespace std;

vector<double> Pairhmm::calculate(const Testcase& testcase) {
  auto results = phmm_float.calculate(testcase);  // calculate all the testcases using the float precision pairhmm implementation
  phmm_double.recalculate(testcase, results, phmm_float.max_original_read_length, phmm_float.max_padded_read_length); // recalculate only the ones that failed with double precision (results = DOUBLE_MIN for failed ones)
  return results;
}

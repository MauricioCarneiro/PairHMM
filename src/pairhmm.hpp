template <typename FAST, typename CORRECT>
std::vector<double> Pairhmm<FAST,CORRECT>::calculate(const Testcase& testcase) {
  auto results = fast.calculate(testcase);  // calculate all the testcases using the float precision pairhmm implementation
  correct.recalculate(testcase, results, fast.max_original_read_length, fast.max_padded_read_length); // recalculate only the ones that failed with double precision (results = DOUBLE_MIN for failed ones)
  return results;
}

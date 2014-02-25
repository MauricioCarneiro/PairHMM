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
    std::vector<double> calculate(const Testcase& testcase);
};

#include "pairhmm.hpp"

#endif

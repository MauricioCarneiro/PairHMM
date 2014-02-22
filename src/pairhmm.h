#ifndef __PAIRHMM__
#define __PAIRHMM__

#include <vector>

#include "pairhmm_impl.h"
#include "testcase.h"

class Pairhmm {
  private:
    PairhmmImpl<float> phmm_float;
    PairhmmImpl<double> phmm_double;

  public: 
    std::vector<double> calculate(const Testcase& testcase);
};

#endif

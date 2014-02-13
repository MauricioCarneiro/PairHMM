#ifndef __TESTCASE__
#define __TESTCASE__

#include <ostream>
#include <vector>

#include "read.h"
#include "haplotype.h"

struct Testcase {
  std::vector<Read<uint8_t>> reads;
  std::vector<Haplotype> haplotypes;

  Testcase() {};
  Testcase(std::vector<Read<uint8_t>>& reads_, std::vector<Haplotype>& haplotypes_) = delete;
  Testcase(std::vector<Read<uint8_t>>&& reads_, std::vector<Haplotype>&& haplotypes_) :
    reads {std::move(reads_)},
    haplotypes {std::move(haplotypes_)}
  {};  

  uint32_t calculate_max_read_length() const;
};

std::ostream& operator<<(std::ostream& out, const Testcase& testcase); 

#endif

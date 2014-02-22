#ifndef __HAPLOTYPE__
#define __HAPLOTYPE__

#include <ostream>
#include <vector>
#include <string>

#include "utils.h"

struct Haplotype {
  size_t original_length;
  std::vector<uint8_t> bases;

  Haplotype(const std::vector<uint8_t> bases_, const size_t original_length_ = 0) : 
    original_length{original_length_},
    bases {bases_}
  {};

  Haplotype(const std::string& bases_, const size_t original_length_ = 0) :
    original_length{original_length_},
    bases {convert_bytes<std::vector<uint8_t>>(bases_)}
  {};
};

std::ostream& operator<<(std::ostream& out, const Haplotype& haplotype); 

#endif

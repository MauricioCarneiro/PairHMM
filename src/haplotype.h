#ifndef __HAPLOTYPE__
#define __HAPLOTYPE__

#include <ostream>
#include <vector>
#include <string>

#include "utils.h"

struct Haplotype {
  std::vector <uint8_t> bases;

  Haplotype(const std::string& bases_) :
    bases {convert_bytes<std::vector<uint8_t>>(bases_)}
  {};

};

std::ostream& operator<<(std::ostream& out, const Haplotype& haplotype); 

#endif

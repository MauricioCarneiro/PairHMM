#include <ostream>

#include "haplotype.h"

std::ostream& operator<<(std::ostream& out, const Haplotype& haplotype) {
  out << convert_bytes<std::string>(haplotype.bases);
  return out;
}

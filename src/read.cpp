#include <ostream>

#include "read.h" 

std::ostream& operator<<(std::ostream& out, const Read& read) {
  out << convert_bytes<std::string>(read.bases) << " ";
  out << convert_bytes<std::string>(read.quals, QUAL_OFFSET) << " ";
  out << convert_bytes<std::string>(read.ins_quals, QUAL_OFFSET) << " ";
  out << convert_bytes<std::string>(read.del_quals, QUAL_OFFSET) << " ";
  out << convert_bytes<std::string>(read.gcp_quals, QUAL_OFFSET) << " ";
  return out;
}

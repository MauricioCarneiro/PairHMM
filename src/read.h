#ifndef __READ__
#define __READ__

#include <ostream>
#include <vector>
#include <string>

#include "utils.h"

struct Read {
  std::vector<uint8_t> bases;
  std::vector<uint8_t> quals;
  std::vector<uint8_t> ins_quals;
  std::vector<uint8_t> del_quals;
  std::vector<uint8_t> gcp_quals;

  Read(const std::string& bases_,
       const std::string& quals_,
       const std::string& ins_quals_,
       const std::string& del_quals_,
       const std::string& gcp_quals_) :
    bases     {convert_bytes<std::vector<uint8_t>>(bases_)},
    quals     {convert_bytes<std::vector<uint8_t>>(quals_)},
    ins_quals {convert_bytes<std::vector<uint8_t>>(ins_quals_, -QUAL_OFFSET)},
    del_quals {convert_bytes<std::vector<uint8_t>>(del_quals_, -QUAL_OFFSET)},
    gcp_quals {convert_bytes<std::vector<uint8_t>>(gcp_quals_, -QUAL_OFFSET)}
  {};

};

std::ostream& operator<<(std::ostream& out, const Read& read);

#endif

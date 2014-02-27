#ifndef __CONSTANTS__
#define __CONSTANTS__

#include <vector> 

#include "read.h"

namespace constants_impl {

template<class PRECISION>
  struct Constants {
    PRECISION mm;
    PRECISION gm;
    PRECISION mx;
    PRECISION xx;
    PRECISION my;
    PRECISION yy;
  };

}

template<class PRECISION>
struct Constants {
  std::vector<constants_impl::Constants<PRECISION>> index;

  Constants(const std::size_t initial_size) :
    index(initial_size) 
  {}

  void resize(const std::size_t new_size) {
    if (index.capacity() < new_size) 
      index.resize(new_size);
  }

  void update(const Read<PRECISION>& read) {
    const auto rows = read.bases.size();
    for (auto i = 1u; i != rows; ++i) {
      index[i].mm = static_cast<PRECISION>(1) - (read.ins_quals[i] + read.del_quals[i]);
      index[i].gm = static_cast<PRECISION>(1) - read.gcp_quals[i];
      index[i].mx = read.ins_quals[i];
      index[i].xx = read.gcp_quals[i];
      index[i].my = i < (rows - 1) ? read.del_quals[i] : static_cast<PRECISION>(1);
      index[i].yy = i < (rows - 1) ? read.gcp_quals[i] : static_cast<PRECISION>(1);
    }
  }

};


#endif

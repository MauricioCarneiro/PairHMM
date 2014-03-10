#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <vector>

#include "read.h"

template<class PRECISION, class ALLOCATOR = std::allocator<PRECISION>>
struct Constants {

  Constants(const std::size_t initial_size): mm(initial_size), gm(initial_size), mx(initial_size),
    xx(initial_size), my(initial_size), yy(initial_size) {}

  void resize(const std::size_t new_size) {
    if (mm.size() < new_size) {
      mm.resize(new_size);
      gm.resize(new_size);
      mx.resize(new_size);
      xx.resize(new_size);
      my.resize(new_size);
      yy.resize(new_size);
    }
  }

  template <class BASE_TYPE>
  void update(const Read<PRECISION, BASE_TYPE>& read) {
    const auto rows = read.bases.size();
    for (auto i = 1u; i != rows; ++i) {
      mm[i] = static_cast<PRECISION>(1) - (read.ins_quals[i] + read.del_quals[i]);
      gm[i] = static_cast<PRECISION>(1) - read.gcp_quals[i];
      mx[i] = read.ins_quals[i];
      xx[i] = read.gcp_quals[i];
      my[i] = i < (rows - 1) ? read.del_quals[i] : static_cast<PRECISION>(1);
      yy[i] = i < (rows - 1) ? read.gcp_quals[i] : static_cast<PRECISION>(1);
    }
  }

  std::vector<PRECISION, ALLOCATOR> mm;
  std::vector<PRECISION, ALLOCATOR> gm;
  std::vector<PRECISION, ALLOCATOR> mx;
  std::vector<PRECISION, ALLOCATOR> xx;
  std::vector<PRECISION, ALLOCATOR> my;
  std::vector<PRECISION, ALLOCATOR> yy;

};


#endif

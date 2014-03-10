#ifndef DIAGONALS_H
#define DIAGONALS_H

#include <vector>
#include <iostream>

template<class PRECISION, class ALLOCATOR = std::allocator<PRECISION>>
struct Diagonals {
  std::vector<PRECISION, ALLOCATOR> m;
  std::vector<PRECISION, ALLOCATOR> mp;
  std::vector<PRECISION, ALLOCATOR> mpp;
  std::vector<PRECISION, ALLOCATOR> x;
  std::vector<PRECISION, ALLOCATOR> xp;
  std::vector<PRECISION, ALLOCATOR> xpp;
  std::vector<PRECISION, ALLOCATOR> y;
  std::vector<PRECISION, ALLOCATOR> yp;
  std::vector<PRECISION, ALLOCATOR> ypp;

  Diagonals() = default;
  Diagonals(const size_t read_length) :
    m   (read_length, static_cast<PRECISION>(0)),
    mp  (read_length, static_cast<PRECISION>(0)),
    mpp (read_length, static_cast<PRECISION>(0)),
    x   (read_length, static_cast<PRECISION>(0)),
    xp  (read_length, static_cast<PRECISION>(0)),
    xpp (read_length, static_cast<PRECISION>(0)),
    y   (read_length, static_cast<PRECISION>(0)),
    yp  (read_length, static_cast<PRECISION>(0)),
    ypp (read_length, static_cast<PRECISION>(0))
  {}

  void resize(const size_t new_size) {
    if (new_size > m.capacity()) { // only grow?
      m.resize(new_size); mp.resize(new_size); mpp.resize(new_size);
      x.resize(new_size); xp.resize(new_size); xpp.resize(new_size);
      y.resize(new_size); yp.resize(new_size); ypp.resize(new_size);
    }
  }

  void rotate() {
    std::swap(mpp, mp);
    std::swap(xpp, xp);
    std::swap(ypp, yp);
    std::swap(mp, m);
    std::swap(xp, x);
    std::swap(yp, y);
  }

  void update(const PRECISION first_row_value) {
    std::fill(m.begin(), m.end(), static_cast<PRECISION>(0)); std::fill(mp.begin(), mp.end(), static_cast<PRECISION>(0)); std::fill(mpp.begin(), mpp.end(), static_cast<PRECISION>(0));
    std::fill(x.begin(), x.end(), static_cast<PRECISION>(0)); std::fill(xp.begin(), xp.end(), static_cast<PRECISION>(0)); std::fill(xpp.begin(), xpp.end(), static_cast<PRECISION>(0));
    std::fill(y.begin(), y.end(), static_cast<PRECISION>(0)); std::fill(yp.begin(), yp.end(), static_cast<PRECISION>(0)); std::fill(ypp.begin(), ypp.end(), static_cast<PRECISION>(0));
    ypp[0] = yp[0] = y[0] = first_row_value;
  }

};

template<class PRECISION>
std::ostream& operator<<(std::ostream& out, const Diagonals<PRECISION>& diagonals);

#endif

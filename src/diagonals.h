#ifndef __DIAGONALS__
#define __DIAGONALS__

#include <vector>
#include <iostream>

template<class precision>
struct Diagonals {
  std::vector<precision> m;
  std::vector<precision> mp;
  std::vector<precision> mpp;
  std::vector<precision> x;
  std::vector<precision> xp;
  std::vector<precision> xpp;
  std::vector<precision> y;
  std::vector<precision> yp;
  std::vector<precision> ypp;

  Diagonals() = default;
  Diagonals(const size_t read_length) :
    m   (read_length, static_cast<precision>(0)),
    mp  (read_length, static_cast<precision>(0)),
    mpp (read_length, static_cast<precision>(0)),
    x   (read_length, static_cast<precision>(0)),
    xp  (read_length, static_cast<precision>(0)),
    xpp (read_length, static_cast<precision>(0)),
    y   (read_length, static_cast<precision>(0)),
    yp  (read_length, static_cast<precision>(0)),
    ypp (read_length, static_cast<precision>(0))
  {}

  void resize(const size_t new_size) {
    if (new_size > m.capacity()) {
      m.resize(new_size); mp.resize(new_size); mpp.resize(new_size);
      x.resize(new_size); xp.resize(new_size); xpp.resize(new_size);
      y.resize(new_size); yp.resize(new_size); ypp.resize(new_size);
    }
  }

  void zero() {
    std::fill(m.begin(), m.end(), static_cast<precision>(0)); std::fill(mp.begin(), mp.end(), static_cast<precision>(0)); std::fill(mpp.begin(), mpp.end(), static_cast<precision>(0));
    std::fill(x.begin(), x.end(), static_cast<precision>(0)); std::fill(xp.begin(), xp.end(), static_cast<precision>(0)); std::fill(xpp.begin(), xpp.end(), static_cast<precision>(0));
    std::fill(y.begin(), y.end(), static_cast<precision>(0)); std::fill(yp.begin(), yp.end(), static_cast<precision>(0)); std::fill(ypp.begin(), ypp.end(), static_cast<precision>(0));
  }

  void swap() {
    std::swap(mpp, mp);
    std::swap(xpp, xp);
    std::swap(ypp, yp);
    std::swap(mp, m); 
    std::swap(xp, x);
    std::swap(yp, y);
  }
    
};

template<class precision>
std::ostream& operator<<(std::ostream& out, const Diagonals<precision>& diagonals); 

#endif

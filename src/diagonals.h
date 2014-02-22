#ifndef __DIAGONALS__
#define __DIAGONALS__

#include <vector>
#include <iostream>

template<class Precision>
struct Diagonals {
  std::vector<Precision> m;
  std::vector<Precision> mp;
  std::vector<Precision> mpp;
  std::vector<Precision> x;
  std::vector<Precision> xp;
  std::vector<Precision> xpp;
  std::vector<Precision> y;
  std::vector<Precision> yp;
  std::vector<Precision> ypp;

  Diagonals() = default;
  Diagonals(const size_t read_length) :
    m   (read_length, static_cast<Precision>(0)),
    mp  (read_length, static_cast<Precision>(0)),
    mpp (read_length, static_cast<Precision>(0)),
    x   (read_length, static_cast<Precision>(0)),
    xp  (read_length, static_cast<Precision>(0)),
    xpp (read_length, static_cast<Precision>(0)),
    y   (read_length, static_cast<Precision>(0)),
    yp  (read_length, static_cast<Precision>(0)),
    ypp (read_length, static_cast<Precision>(0))
  {}

  void resize(const size_t new_size) {
    if (new_size > m.capacity()) {
      m.resize(new_size); mp.resize(new_size); mpp.resize(new_size);
      x.resize(new_size); xp.resize(new_size); xpp.resize(new_size);
      y.resize(new_size); yp.resize(new_size); ypp.resize(new_size);
    }
  }

  void zero() {
    std::fill(m.begin(), m.end(), static_cast<Precision>(0)); std::fill(mp.begin(), mp.end(), static_cast<Precision>(0)); std::fill(mpp.begin(), mpp.end(), static_cast<Precision>(0));
    std::fill(x.begin(), x.end(), static_cast<Precision>(0)); std::fill(xp.begin(), xp.end(), static_cast<Precision>(0)); std::fill(xpp.begin(), xpp.end(), static_cast<Precision>(0));
    std::fill(y.begin(), y.end(), static_cast<Precision>(0)); std::fill(yp.begin(), yp.end(), static_cast<Precision>(0)); std::fill(ypp.begin(), ypp.end(), static_cast<Precision>(0));
  }

  void rotate() {
    std::swap(mpp, mp);
    std::swap(xpp, xp);
    std::swap(ypp, yp);
    std::swap(mp, m); 
    std::swap(xp, x);
    std::swap(yp, y);
  }
    
};

template<class Precision>
std::ostream& operator<<(std::ostream& out, const Diagonals<Precision>& diagonals); 

#endif

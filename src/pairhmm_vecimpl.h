#ifndef PAIRHMM_VECIMPL_H
#define PAIRHMM_VECIMPL_H

#include "pairhmm_impl.h"

template <class PRECISION, class DIAGONALS, class CONSTANTS, int VECSIZE = 1>
class PairhmmVecImpl: public PairhmmImpl<PRECISION, DIAGONALS, CONSTANTS, VECSIZE> {
 using Base = PairhmmImpl<PRECISION,DIAGONALS,CONSTANTS,VECSIZE>;
 public:
  PairhmmVecImpl(const size_t initial_size = Base::INITIAL_SIZE) : Base {initial_size} {}
  virtual ~PairhmmVecImpl() { }
 protected:
  virtual double do_compute_full_prob(const Read<PRECISION,PRECISION>& read, const Haplotype<PRECISION>& haplotype) override {
    const auto hl = haplotype.original_length;  // haplotype original length (unpadded)
    const auto rl = read.original_length;       // read original length (unpadded)
    const auto rows = rl + read.left_padding;  // number of rows in the diagonals (padded read length)
    const auto mrl = this->max_original_read_length();  // alias for max original read length for readability in the code below (max read length in the testcase)
    const auto fd = mrl - rl;                   // first diagonal to compute (saves compute of all-0 diagonals when read is shorter than the padding - which will be the maximum read length in the testcase)
    auto result = 0.l;                          // result accumulator
    auto &diags = this->m_diagonals;
    auto &consts = this->m_constants;
    auto N = PRECISION('N');
    const auto hap_last = mrl+hl-1;
    for (auto d = fd; d != hap_last; ++d) { // d for diagonal
      for (auto r = 1u; r < rows; r += VECSIZE) {       // r for row
        for (auto v = 0; v < VECSIZE; v++) {
          const auto read_base = read.bases[r+v];
          const auto hap_base = haplotype.bases[hap_last+r+v-d];
          const auto base_qual = read.base_quals[r+v];
          const auto prior = ((read_base == hap_base) || (read_base == N) || (hap_base == N)) ?  PRECISION(1) - base_qual : base_qual;
          diags.m[r+v] = prior * ((diags.mpp[r+v-1] * consts.mm[r+v]) + (consts.gm[r+v] * (diags.xpp[r+v-1] + diags.ypp[r+v-1])));
          diags.x[r+v] = diags.mp[r+v-1] * consts.mx[r+v] + diags.xp[r+v-1] * consts.xx[r+v];
          diags.y[r+v] = diags.mp[r+v] * consts.my[r+v] + diags.yp[r+v] * consts.yy[r+v];
        }
      }
      result += diags.m[rows-1] + diags.x[rows-1];
      diags.rotate();
    }
    return result < this->MIN_ACCEPTED ?
      this->FAILED_RUN_RESULT :                       // if we underflowed return failed constant to rerun with higher precision if desired
      log10(static_cast<double>(result)) - log10(static_cast<double>(this->INITIAL_CONSTANT));
  }
};

#endif

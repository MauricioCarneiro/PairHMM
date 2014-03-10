#ifndef PAIRHMM_SSEIMPL_H
#define PAIRHMM_SSEIMPL_H

#include <xmmintrin.h>
#include "aligned_allocator.h"
#include "pairhmm_impl.h"

class PairhmmSSEFloatImpl: public PairhmmImpl<float, Diagonals3<float, Aligned_allocator<float, 16, 4>>, Constants<float, Aligned_allocator<float, 16, 4>>, 4> {
 using Base = PairhmmImpl<float, Diagonals3<float, Aligned_allocator<float, 16, 4>>, Constants<float, Aligned_allocator<float, 16, 4>>, 4>;
public:
  PairhmmSSEFloatImpl(const size_t initial_size = Base::INITIAL_SIZE) : Base {initial_size} { }
  virtual ~PairhmmSSEFloatImpl() { }
protected:
  virtual double do_compute_full_prob(const Read<float,float>& read, const Haplotype<float>& haplotype) override {
    const auto hl = haplotype.original_length;  // haplotype original length (unpadded)
    const auto rl = read.original_length;       // read original length (unpadded)
    const auto rows = rl + read.left_padding;  // number of rows in the diagonals (padded read length)
    const auto mrl = this->max_original_read_length();  // alias for max original read length for readability in the code below (max read length in the testcase)
    const auto fd = mrl - rl;                   // first diagonal to compute (saves compute of all-0 diagonals when read is shorter than the padding - which will be the maximum read length in the testcase)
    auto result = 0.l;                          // result accumulator
    auto &diags = this->m_diagonals;
    auto &consts = this->m_constants;
    const auto hap_last = mrl+hl-1;
    const __m128i N = _mm_castps_si128(_mm_set1_ps('N'));
    const __m128 one = _mm_set1_ps(1.f);
    for (auto d = fd; d != hap_last; ++d) { // d for diagonal
      for (auto r = 1u; r < rows; r += 4) {       // r for row
        __m128i read_base = _mm_castps_si128(_mm_loadu_ps(&read.bases[r]));
        __m128i hap_base = _mm_castps_si128(_mm_loadu_ps(&haplotype.bases[hap_last+r-d]));
        __m128 base_qual = _mm_loadu_ps(&read.base_quals[r]);
        __m128 one_minus_base_qual = _mm_sub_ps(one, base_qual);
        __m128i cmp = _mm_or_si128(_mm_cmpeq_epi32(read_base, hap_base),
                        _mm_or_si128(_mm_cmpeq_epi32(read_base, N), _mm_cmpeq_epi32(hap_base, N)));

        __m128 prior = _mm_castsi128_ps(
            _mm_or_si128(
                _mm_and_si128(cmp, _mm_castps_si128(one_minus_base_qual)),
                _mm_andnot_si128(cmp, _mm_castps_si128(base_qual))));

        _mm_store_ps(&diags.m[r], _mm_mul_ps(prior,
            _mm_add_ps(_mm_mul_ps(_mm_loadu_ps(&diags.mpp[r-1]), _mm_loadu_ps(&consts.mm[r])),
                _mm_mul_ps(_mm_loadu_ps(&consts.gm[r]),
                    _mm_add_ps(_mm_loadu_ps(&diags.xpp[r-1]), _mm_loadu_ps(&diags.ypp[r-1]))))));

        _mm_store_ps(&diags.x[r], _mm_add_ps(
            _mm_mul_ps(_mm_loadu_ps(&diags.mp[r-1]), _mm_loadu_ps(&consts.mx[r])),
            _mm_mul_ps(_mm_loadu_ps(&diags.xp[r-1]), _mm_loadu_ps(&consts.xx[r]))));

        _mm_store_ps(&diags.y[r], _mm_add_ps(
            _mm_mul_ps(_mm_loadu_ps(&diags.mp[r]), _mm_loadu_ps(&consts.my[r])),
            _mm_mul_ps(_mm_loadu_ps(&diags.yp[r]), _mm_loadu_ps(&consts.yy[r]))));
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

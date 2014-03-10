#ifndef PAIRHMM_AVXIMPL_H
#define PAIRHMM_AVXIMPL_H

#include <immintrin.h>
#include "aligned_allocator.h"
#include "pairhmm_impl.h"

class PairhmmAVXFloatImpl: public PairhmmImpl<float, Diagonals<float, Aligned_allocator<float, 32, 4>>, Constants<float, Aligned_allocator<float, 32, 4>>, 8> {
 using Base = PairhmmImpl<float, Diagonals<float, Aligned_allocator<float, 32, 4>>, Constants<float, Aligned_allocator<float, 32, 4>>, 8>;
public:
  PairhmmAVXFloatImpl(const size_t initial_size = Base::INITIAL_SIZE) : Base {initial_size} {
    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
    _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
  }
  virtual ~PairhmmAVXFloatImpl() { }
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
    const __m256 N = _mm256_set1_ps('N');
    const __m256 one = _mm256_set1_ps(1.f);
    for (auto d = fd; d != hap_last; ++d) { // d for diagonal
      for (auto r = 1u; r < rows; r += 8) {       // r for row
        __m256 read_base = _mm256_loadu_ps(&read.bases[r]);
        __m256 hap_base = _mm256_loadu_ps(&haplotype.bases[hap_last+r-d]);
        __m256 base_qual = _mm256_loadu_ps(&read.base_quals[r]);
        __m256 one_minus_base_qual = _mm256_sub_ps(one, base_qual);
        __m256 cmp = _mm256_or_ps(_mm256_cmp_ps(read_base, hap_base, _CMP_EQ_UQ),
          _mm256_or_ps(_mm256_cmp_ps(read_base, N, _CMP_EQ_UQ), _mm256_cmp_ps(hap_base, N, _CMP_EQ_UQ)));
        __m256 prior = _mm256_blendv_ps(base_qual, one_minus_base_qual, cmp);
        _mm256_store_ps(&diags.m[r], _mm256_mul_ps(prior,
            _mm256_add_ps(_mm256_mul_ps(_mm256_loadu_ps(&diags.mpp[r-1]), _mm256_loadu_ps(&consts.mm[r])),
                _mm256_mul_ps(_mm256_loadu_ps(&consts.gm[r]),
                    _mm256_add_ps(_mm256_loadu_ps(&diags.xpp[r-1]), _mm256_loadu_ps(&diags.ypp[r-1]))))));
        _mm256_store_ps(&diags.x[r], _mm256_add_ps(
            _mm256_mul_ps(_mm256_loadu_ps(&diags.mp[r-1]), _mm256_loadu_ps(&consts.mx[r])),
            _mm256_mul_ps(_mm256_loadu_ps(&diags.xp[r-1]), _mm256_loadu_ps(&consts.xx[r]))));
        _mm256_store_ps(&diags.y[r], _mm256_add_ps(
            _mm256_mul_ps(_mm256_loadu_ps(&diags.mp[r]), _mm256_loadu_ps(&consts.my[r])),
            _mm256_mul_ps(_mm256_loadu_ps(&diags.yp[r]), _mm256_loadu_ps(&consts.yy[r]))));
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

#ifndef PAIRHMM_AVXIMPL_H
#define PAIRHMM_AVXIMPL_H

#include <immintrin.h>
#include "aligned_allocator.h"
#include "pairhmm_impl.h"

class PairhmmAVXFloatImpl: public PairhmmImpl<float, Diagonals3<float, Aligned_allocator<float, 32, 4>>, Constants<float, Aligned_allocator<float, 32, 4>>, 8> {
 using Base = PairhmmImpl<float, Diagonals3<float, Aligned_allocator<float, 32, 4>>, Constants<float, Aligned_allocator<float, 32, 4>>, 8>;
public:
  PairhmmAVXFloatImpl(const size_t initial_size = Base::INITIAL_SIZE) : Base {initial_size} { }
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
            _mm256_add_ps(_mm256_mul_ps(_mm256_loadu_ps(&diags.mpp[r-1]), _mm256_load_ps(&consts.mm[r])),
                _mm256_mul_ps(_mm256_load_ps(&consts.gm[r]),
                    _mm256_add_ps(_mm256_loadu_ps(&diags.xpp[r-1]), _mm256_loadu_ps(&diags.ypp[r-1]))))));
        _mm256_store_ps(&diags.x[r], _mm256_add_ps(
            _mm256_mul_ps(_mm256_loadu_ps(&diags.mp[r-1]), _mm256_load_ps(&consts.mx[r])),
            _mm256_mul_ps(_mm256_loadu_ps(&diags.xp[r-1]), _mm256_load_ps(&consts.xx[r]))));
        _mm256_store_ps(&diags.y[r], _mm256_add_ps(
            _mm256_mul_ps(_mm256_load_ps(&diags.mp[r]), _mm256_load_ps(&consts.my[r])),
            _mm256_mul_ps(_mm256_load_ps(&diags.yp[r]), _mm256_load_ps(&consts.yy[r]))));
      }
      result += diags.m[rows-1] + diags.x[rows-1];
      diags.rotate();
    }
    return result < this->MIN_ACCEPTED ?
      this->FAILED_RUN_RESULT :                       // if we underflowed return failed constant to rerun with higher precision if desired
      log10(static_cast<double>(result)) - log10(static_cast<double>(this->INITIAL_CONSTANT));
  }
};

class PairhmmAVXDoubleImpl: public PairhmmImpl<double, Diagonals3<double, Aligned_allocator<double, 32, 8>>, Constants<double, Aligned_allocator<double, 32, 8>>, 4> {
 using Base = PairhmmImpl<double, Diagonals3<double, Aligned_allocator<double, 32, 8>>, Constants<double, Aligned_allocator<double, 32, 8>>, 4>;
public:
  PairhmmAVXDoubleImpl(const size_t initial_size = Base::INITIAL_SIZE) : Base {initial_size} { }
  virtual ~PairhmmAVXDoubleImpl() { }
protected:
  virtual double do_compute_full_prob(const Read<double,double>& read, const Haplotype<double>& haplotype) override {
    const auto hl = haplotype.original_length;  // haplotype original length (unpadded)
    const auto rl = read.original_length;       // read original length (unpadded)
    const auto rows = rl + read.left_padding;  // number of rows in the diagonals (padded read length)
    const auto mrl = this->max_original_read_length();  // alias for max original read length for readability in the code below (max read length in the testcase)
    const auto fd = mrl - rl;                   // first diagonal to compute (saves compute of all-0 diagonals when read is shorter than the padding - which will be the maximum read length in the testcase)
    auto result = 0.l;                          // result accumulator
    auto &diags = this->m_diagonals;
    auto &consts = this->m_constants;
    const auto hap_last = mrl+hl-1;
    const __m256d N = _mm256_set1_pd('N');
    const __m256d one = _mm256_set1_pd(1.);
    for (auto d = fd; d != hap_last; ++d) { // d for diagonal
      for (auto r = 1u; r < rows; r += 4) {       // r for row
        __m256d read_base = _mm256_loadu_pd(&read.bases[r]);
        __m256d hap_base = _mm256_loadu_pd(&haplotype.bases[hap_last+r-d]);
        __m256d base_qual = _mm256_loadu_pd(&read.base_quals[r]);
        __m256d one_minus_base_qual = _mm256_sub_pd(one, base_qual);
        __m256d cmp = _mm256_or_pd(_mm256_cmp_pd(read_base, hap_base, _CMP_EQ_UQ),
          _mm256_or_pd(_mm256_cmp_pd(read_base, N, _CMP_EQ_UQ), _mm256_cmp_pd(hap_base, N, _CMP_EQ_UQ)));
        __m256d prior = _mm256_blendv_pd(base_qual, one_minus_base_qual, cmp);
        _mm256_store_pd(&diags.m[r], _mm256_mul_pd(prior,
            _mm256_add_pd(_mm256_mul_pd(_mm256_loadu_pd(&diags.mpp[r-1]), _mm256_load_pd(&consts.mm[r])),
                _mm256_mul_pd(_mm256_load_pd(&consts.gm[r]),
                    _mm256_add_pd(_mm256_loadu_pd(&diags.xpp[r-1]), _mm256_loadu_pd(&diags.ypp[r-1]))))));
        _mm256_store_pd(&diags.x[r], _mm256_add_pd(
            _mm256_mul_pd(_mm256_loadu_pd(&diags.mp[r-1]), _mm256_load_pd(&consts.mx[r])),
            _mm256_mul_pd(_mm256_loadu_pd(&diags.xp[r-1]), _mm256_load_pd(&consts.xx[r]))));
        _mm256_store_pd(&diags.y[r], _mm256_add_pd(
            _mm256_mul_pd(_mm256_load_pd(&diags.mp[r]), _mm256_load_pd(&consts.my[r])),
            _mm256_mul_pd(_mm256_load_pd(&diags.yp[r]), _mm256_load_pd(&consts.yy[r]))));
      }
      result += diags.m[rows-1] + diags.x[rows-1];
      diags.rotate();
    }
    return result < this->MIN_ACCEPTED ?
      this->FAILED_RUN_RESULT :                       // if we underflowed return failed constant to rerun with higher precision if desired
      log10(static_cast<double>(result)) - log10(static_cast<double>(this->INITIAL_CONSTANT));
  }
};

class PairhmmAVXFloat2DiagsImpl: public PairhmmImpl<float, Diagonals2<float, Aligned_allocator<float, 32, 4>>, Constants<float, Aligned_allocator<float, 32, 4>>, 8> {
 using Base = PairhmmImpl<float, Diagonals2<float, Aligned_allocator<float, 32, 4>>, Constants<float, Aligned_allocator<float, 32, 4>>, 8>;
public:
  PairhmmAVXFloat2DiagsImpl(const size_t initial_size = Base::INITIAL_SIZE) : Base {initial_size} { }
  virtual ~PairhmmAVXFloat2DiagsImpl() { }
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
      __m256 xpp = _mm256_loadu_ps(&diags.x[0]); // Warning: magic follows.
      __m256 ypp = _mm256_loadu_ps(&diags.y[0]);
      __m256 mpp = _mm256_loadu_ps(&diags.m[0]);
      for (auto r = 1u; r < rows; r += 8) {       // r for row
        __m256 read_base = _mm256_loadu_ps(&read.bases[r]); // maybe we should align read.bases and read.base_quals??
        __m256 hap_base = _mm256_loadu_ps(&haplotype.bases[hap_last+r-d]);
        __m256 base_qual = _mm256_loadu_ps(&read.base_quals[r]);
        __m256 one_minus_base_qual = _mm256_sub_ps(one, base_qual);
        __m256 cmp = _mm256_or_ps(_mm256_cmp_ps(read_base, hap_base, _CMP_EQ_UQ),
            _mm256_or_ps(_mm256_cmp_ps(read_base, N, _CMP_EQ_UQ), _mm256_cmp_ps(hap_base, N, _CMP_EQ_UQ)));
        __m256 prior = _mm256_blendv_ps(base_qual, one_minus_base_qual, cmp);
        __m256 mpp_tmp = _mm256_loadu_ps(&diags.m[r+8-1]);
        _mm256_store_ps(&diags.m[r], _mm256_mul_ps(prior,
            _mm256_add_ps(_mm256_mul_ps(mpp, _mm256_load_ps(&consts.mm[r])),
                _mm256_mul_ps(_mm256_load_ps(&consts.gm[r]), _mm256_add_ps(xpp, ypp)))));
        mpp = mpp_tmp;
        xpp = _mm256_loadu_ps(&diags.x[r+8-1]);
        _mm256_store_ps(&diags.x[r], _mm256_add_ps(
            _mm256_mul_ps(_mm256_loadu_ps(&diags.mp[r-1]), _mm256_load_ps(&consts.mx[r])),
            _mm256_mul_ps(_mm256_loadu_ps(&diags.xp[r-1]), _mm256_load_ps(&consts.xx[r]))));
        ypp = _mm256_loadu_ps(&diags.y[r+8-1]);
        _mm256_store_ps(&diags.y[r], _mm256_add_ps(
            _mm256_mul_ps(_mm256_load_ps(&diags.mp[r]), _mm256_load_ps(&consts.my[r])),
            _mm256_mul_ps(_mm256_load_ps(&diags.yp[r]), _mm256_load_ps(&consts.yy[r]))));
      }
      result += diags.m[rows-1] + diags.x[rows-1];
      diags.rotate();
    }
    return result < this->MIN_ACCEPTED ?
      this->FAILED_RUN_RESULT :                       // if we underflowed return failed constant to rerun with higher precision if desired
      log10(static_cast<double>(result)) - log10(static_cast<double>(this->INITIAL_CONSTANT));
  }
};


class PairhmmAVXDouble2DiagsImpl: public PairhmmImpl<double, Diagonals2<double, Aligned_allocator<double, 32, 8>>, Constants<double, Aligned_allocator<double, 32, 8>>, 4> {
 using Base = PairhmmImpl<double, Diagonals2<double, Aligned_allocator<double, 32, 8>>, Constants<double, Aligned_allocator<double, 32, 8>>, 4>;
public:
  PairhmmAVXDouble2DiagsImpl(const size_t initial_size = Base::INITIAL_SIZE) : Base {initial_size} { }
  virtual ~PairhmmAVXDouble2DiagsImpl() { }
protected:
  virtual double do_compute_full_prob(const Read<double,double>& read, const Haplotype<double>& haplotype) override {
    const auto hl = haplotype.original_length;  // haplotype original length (unpadded)
    const auto rl = read.original_length;       // read original length (unpadded)
    const auto rows = rl + read.left_padding;  // number of rows in the diagonals (padded read length)
    const auto mrl = this->max_original_read_length();  // alias for max original read length for readability in the code below (max read length in the testcase)
    const auto fd = mrl - rl;                   // first diagonal to compute (saves compute of all-0 diagonals when read is shorter than the padding - which will be the maximum read length in the testcase)
    auto result = 0.l;                          // result accumulator
    auto &diags = this->m_diagonals;
    auto &consts = this->m_constants;
    const auto hap_last = mrl+hl-1;
    const __m256d N = _mm256_set1_pd('N');
    const __m256d one = _mm256_set1_pd(1.);
    for (auto d = fd; d != hap_last; ++d) { // d for diagonal
      __m256d xpp = _mm256_loadu_pd(&diags.x[0]); // Warning: magic follows.
      __m256d ypp = _mm256_loadu_pd(&diags.y[0]);
      __m256d mpp = _mm256_loadu_pd(&diags.m[0]);
      for (auto r = 1u; r < rows; r += 4) {       // r for row
        __m256d read_base = _mm256_loadu_pd(&read.bases[r]);
        __m256d hap_base = _mm256_loadu_pd(&haplotype.bases[hap_last+r-d]);
        __m256d base_qual = _mm256_loadu_pd(&read.base_quals[r]);
        __m256d one_minus_base_qual = _mm256_sub_pd(one, base_qual);
        __m256d cmp = _mm256_or_pd(_mm256_cmp_pd(read_base, hap_base, _CMP_EQ_UQ),
            _mm256_or_pd(_mm256_cmp_pd(read_base, N, _CMP_EQ_UQ), _mm256_cmp_pd(hap_base, N, _CMP_EQ_UQ)));
        __m256d prior = _mm256_blendv_pd(base_qual, one_minus_base_qual, cmp);
        __m256d mpp_tmp = _mm256_loadu_pd(&diags.m[r+4-1]);
        _mm256_store_pd(&diags.m[r], _mm256_mul_pd(prior,
            _mm256_add_pd(_mm256_mul_pd(mpp, _mm256_load_pd(&consts.mm[r])),
                _mm256_mul_pd(_mm256_load_pd(&consts.gm[r]), _mm256_add_pd(xpp, ypp)))));
        mpp = mpp_tmp;
        xpp = _mm256_loadu_pd(&diags.x[r+4-1]);
        _mm256_store_pd(&diags.x[r], _mm256_add_pd(
            _mm256_mul_pd(_mm256_loadu_pd(&diags.mp[r-1]), _mm256_load_pd(&consts.mx[r])),
            _mm256_mul_pd(_mm256_loadu_pd(&diags.xp[r-1]), _mm256_load_pd(&consts.xx[r]))));
        ypp = _mm256_loadu_pd(&diags.y[r+4-1]);
        _mm256_store_pd(&diags.y[r], _mm256_add_pd(
            _mm256_mul_pd(_mm256_load_pd(&diags.mp[r]), _mm256_load_pd(&consts.my[r])),
            _mm256_mul_pd(_mm256_load_pd(&diags.yp[r]), _mm256_load_pd(&consts.yy[r]))));
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

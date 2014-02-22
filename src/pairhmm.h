#ifndef __PAIRHMM__
#define __PAIRHMM__
#include <vector>
#include <cmath>
#include <algorithm>
#include <ostream>

#include "testcase.h"
#include "constants.h"
#include "diagonals.h"

template <class precision>
class Pairhmm {
 private:
  std::vector<Constants<precision>> constants;
  Diagonals<precision> diagonals;
  std::vector<precision> ph2pr;

  static constexpr auto MAX_VAL = 128;
  static constexpr auto INITIAL_CONSTANT = 1e32f;
  static constexpr auto INITIAL_SIZE = 25;

  template <class T>
  void pad(std::vector<T>& v, const uint32_t leading_pad = 1, const uint32_t trailing_pad = 0) const {
    auto total_padding = leading_pad + trailing_pad;
    v.resize(v.size() + total_padding);
    copy(v.begin(), v.end() - total_padding, v.begin() + leading_pad);
  }

  std::vector<Read<float>> pad_reads(std::vector<Read<uint8_t>>& reads) const {
    auto padded_reads = std::vector<Read<float>>{};
    for (auto& read : reads) 
      padded_reads.push_back(pad_read(read));
    return padded_reads;
  }

  std::vector<Haplotype> pad_haplotypes(std::vector<Haplotype>& haplotypes, const size_t max_read_length) const {
    auto padded_haplotypes = std::vector<Haplotype>{};
    for (auto& haplotype : haplotypes) 
      padded_haplotypes.push_back(pad_haplotype(haplotype, max_read_length));
    return padded_haplotypes;
  }

  std::vector<float> real(const std::vector<uint8_t>& v) const {
    auto result = std::vector<float>{};
    for (const auto& q : v)
      result.push_back(ph2pr[q]);
    return result;
  }

  const Haplotype& pad_haplotype(Haplotype& haplotype, const uint32_t read_length) const {
    haplotype.original_length = haplotype.bases.size();  
    std::reverse(haplotype.bases.begin(), haplotype.bases.end());
    pad(haplotype.bases, read_length + 1, read_length);
    return haplotype;
  }

  Read<float> pad_read(const Read<uint8_t>& read) const {
    const auto rl = read.bases.size();
    auto padded_read = Read<float>{read.bases,real(read.base_quals),real(read.ins_quals),real(read.del_quals),real(read.gcp_quals),rl};
    pad(padded_read.bases);
    pad(padded_read.base_quals);
    pad(padded_read.ins_quals);
    pad(padded_read.del_quals);
    pad(padded_read.gcp_quals);
    return padded_read;
  }

  void resize_constants(const size_t new_size) {
    if (constants.capacity() < new_size) 
      constants.resize(new_size);
  }

  template <class T>
  size_t calculate_max_read_length(const std::vector<Read<T>>& reads) const {
    auto max_read_length = 0u;
    for (const auto& read : reads)
      max_read_length = max_read_length >= read.bases.size() ? max_read_length : read.bases.size();
    return max_read_length;
  }

  void update_constants(const Read<float>& read) {
    const auto rows = read.bases.size();
    for (auto i = 1u; i != rows; ++i) {
      constants[i].mm = 1.f - (read.ins_quals[i] + read.del_quals[i]);
      constants[i].gm = 1.f - read.gcp_quals[i];
      constants[i].mx = read.ins_quals[i];
      constants[i].xx = read.gcp_quals[i];
      constants[i].my = i < (rows - 1) ? read.del_quals[i] : 1.f;
      constants[i].yy = i < (rows - 1) ? read.gcp_quals[i] : 1.f;
    }
  }

  void update_diagonals(const Haplotype& hap) {
    diagonals.zero();
    const auto first_row_value = INITIAL_CONSTANT / hap.original_length;
    diagonals.ypp[0] = diagonals.yp[0] = diagonals.y[0] = first_row_value;
  }

  inline float calculate_prior(const char read_base, const char hap_base, const float base_qual) const {
    return  ((read_base == hap_base) || (read_base == 'N') || (hap_base == 'N')) ? 
      1.f - base_qual : 
      base_qual;
  }

  double compute_full_prob(const Read<float>& read, const Haplotype& haplotype) {
    const auto rows = read.bases.size();
    const auto rl = read.original_length;
    const auto hl = haplotype.original_length;
    auto result = 0.l;
    for (auto d = 0u; d != rl + hl - 1; ++d) {  // d for diagonal
      for (auto r = 1u; r != rows; ++r) {  // r for row
        auto hap_idx = rl+hl-d+r-1;
        auto read_base = read.bases[r];
        auto hap_base = haplotype.bases[hap_idx];
        auto prior = calculate_prior(read_base, hap_base, read.base_quals[r]);
        diagonals.m[r] = prior * ((diagonals.mpp[r-1] * constants[r].mm) + (constants[r].gm * (diagonals.xpp[r-1] + diagonals.ypp[r-1])));
        diagonals.x[r] = diagonals.mp[r-1] * constants[r].mx + diagonals.xp[r-1] * constants[r].xx; 
        diagonals.y[r] = diagonals.mp[r] * constants[r].my + diagonals.yp[r] * constants[r].yy; 
      }
      result += diagonals.m[rows-1] + diagonals.x[rows-1];
      diagonals.swap();
    }
    return log10(static_cast<double>(result)) - log10(static_cast<double>(INITIAL_CONSTANT));
  }

  std::vector<double> calculate (const std::vector<Read<float>>& padded_reads, const std::vector<Haplotype>& padded_haplotypes) {
    auto results = std::vector<double>{};
    for (const auto& read : padded_reads) {
      update_constants(read);
      for (const auto& hap : padded_haplotypes) {
        update_diagonals(hap);
        results.push_back(compute_full_prob(read, hap));
      }
    }
    return results;
  }

 public:

  Pairhmm(const size_t initial_size = INITIAL_SIZE) :
    constants {initial_size},
    diagonals {initial_size},
    ph2pr{}
  {
    for (auto i=0.f; i!=MAX_VAL; ++i) 
      ph2pr.push_back(pow(static_cast<precision>(10.0), (-i) / static_cast<precision>(10.0))); 
  }

  std::vector<double> calculate (Testcase& testcase) {
    const auto max_original_read_length = calculate_max_read_length(testcase.reads);
    const auto padded_haplotypes = pad_haplotypes(testcase.haplotypes, max_original_read_length);
    const auto padded_reads = pad_reads(testcase.reads);
    const auto max_padded_read_length = calculate_max_read_length(padded_reads);
    resize_constants(max_padded_read_length);
    diagonals.resize(max_padded_read_length);
    return calculate(padded_reads, padded_haplotypes);
  }

};

#endif

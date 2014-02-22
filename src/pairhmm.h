#ifndef __PAIRHMM__
#define __PAIRHMM__
#include <vector>
#include <cmath>
#include <algorithm>
#include <ostream>

#include "testcase.h"
#include "constants.h"
#include "diagonals.h"

template <class Precision>
class Pairhmm {
 private:
  std::vector<Constants<Precision>> constants;
  Diagonals<Precision> diagonals;
  std::vector<Precision> ph2pr;

  static constexpr auto MAX_PH2PR_INDEX = 128;
  static constexpr auto INITIAL_CONSTANT = 1e32f;
  static constexpr auto INITIAL_SIZE = 250;

  template<class T>
  void pad (T& v, size_t padding) const {
    while(padding-- > 0)
      v.push_back(0);
  }

  Haplotype pad_and_reverse_haplotype(const Haplotype& haplotype, const size_t left_padding, const size_t right_padding) const {
    const auto total_padding = left_padding + haplotype.bases.size() + right_padding;
    auto p = Haplotype{};
    p.original_length = haplotype.bases.size();  
    p.bases.reserve(total_padding);
    pad(p.bases, left_padding);
    for (auto i = haplotype.bases.rbegin(); i != haplotype.bases.rend(); ++i)
      p.bases.push_back(*i);
    pad(p.bases, right_padding);
    return p;
  }

  const Haplotype pad_haplotype(const Haplotype& haplotype, const uint32_t read_length) const {
    const auto left_padding = read_length + 1;
    const auto right_padding = read_length;
    const auto padded_haplotype = pad_and_reverse_haplotype(haplotype, left_padding, right_padding);
    return padded_haplotype;
  }

  template<class Result_Type>
  std::vector<Result_Type> pad_and_convert_qual(const std::vector<uint8_t>& qual, const size_t left_padding, const size_t right_padding = 0, const bool convert = true) const {
    const auto total_padding = left_padding + qual.size() + right_padding;
    auto p = std::vector<Result_Type>{};
    p.reserve(total_padding);
    pad(p, left_padding);
    for (const auto q : qual)
      p.push_back(convert ? ph2pr[q] : q);
    pad(p, right_padding);
    return p;
  }

  Read<Precision> pad_read(const Read<uint8_t>& read) const {
    const auto left_padding = 1;
    auto padded_read = Read<Precision>{};
    padded_read.original_length = read.bases.size();
    padded_read.bases      = pad_and_convert_qual<uint8_t>(read.bases, left_padding, 0, false);
    padded_read.base_quals = pad_and_convert_qual<Precision>(read.base_quals, left_padding);
    padded_read.ins_quals  = pad_and_convert_qual<Precision>(read.ins_quals, left_padding);
    padded_read.del_quals  = pad_and_convert_qual<Precision>(read.del_quals, left_padding);
    padded_read.gcp_quals  = pad_and_convert_qual<Precision>(read.gcp_quals, left_padding);
    return padded_read;
  }

  std::vector<Read<Precision>> pad_reads(std::vector<Read<uint8_t>>& reads) const {
    auto padded_reads = std::vector<Read<Precision>>{};
    for (auto& read : reads) 
      padded_reads.push_back(pad_read(read));
    return padded_reads;
  }

  std::vector<Haplotype> pad_haplotypes(const std::vector<Haplotype>& haplotypes, const size_t max_read_length) const {
    auto padded_haplotypes = std::vector<Haplotype>{};
    padded_haplotypes.reserve(haplotypes.size());
    for (const auto& haplotype : haplotypes) 
      padded_haplotypes.push_back(pad_haplotype(haplotype, max_read_length));
    return padded_haplotypes;
  }

  template <class T>
  size_t calculate_max_read_length(const std::vector<Read<T>>& reads) const {
    auto max_read_length = 0u;
    for (const auto& read : reads)
      max_read_length = max_read_length >= read.bases.size() ? max_read_length : read.bases.size();
    return max_read_length;
  }

  void resize_constants(const size_t new_size) {
    if (constants.capacity() < new_size) 
      constants.resize(new_size);
  }

  void update_constants(const Read<Precision>& read) {
    const auto rows = read.bases.size();
    for (auto i = 1u; i != rows; ++i) {
      constants[i].mm = static_cast<Precision>(1) - (read.ins_quals[i] + read.del_quals[i]);
      constants[i].gm = static_cast<Precision>(1) - read.gcp_quals[i];
      constants[i].mx = read.ins_quals[i];
      constants[i].xx = read.gcp_quals[i];
      constants[i].my = i < (rows - 1) ? read.del_quals[i] : static_cast<Precision>(1);
      constants[i].yy = i < (rows - 1) ? read.gcp_quals[i] : static_cast<Precision>(1);
    }
  }

  void update_diagonals(const Haplotype& hap) {
    diagonals.zero();
    const auto first_row_value = INITIAL_CONSTANT / hap.original_length;
    diagonals.ypp[0] = diagonals.yp[0] = diagonals.y[0] = first_row_value;
  }

  inline Precision calculate_prior(const char read_base, const char hap_base, const Precision base_qual) const {
    return  ((read_base == hap_base) || (read_base == 'N') || (hap_base == 'N')) ? 
      static_cast<Precision>(1) - base_qual : 
      base_qual;
  }

  double compute_full_prob(const Read<Precision>& read, const Haplotype& haplotype) {
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
      diagonals.rotate();
    }
    return log10(static_cast<double>(result)) - log10(static_cast<double>(INITIAL_CONSTANT));
  }

  std::vector<double> calculate (const std::vector<Read<Precision>>& padded_reads, const std::vector<Haplotype>& padded_haplotypes) {
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
    for (auto i=static_cast<Precision>(0); i!=MAX_PH2PR_INDEX; ++i) 
      ph2pr.push_back(pow(static_cast<Precision>(10.0), (-i) / static_cast<Precision>(10.0))); 
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

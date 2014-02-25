#ifndef __PAIRHMM_IMPL__
#define __PAIRHMM_IMPL__

#include <vector>
#include <cmath>
#include <algorithm>
#include <ostream>
#include <unordered_map>
#include <cassert>

#include "testcase.h"
#include "constants.h"
#include "diagonals.h"

#define P(x) x 

namespace constants_with_precision {
  template <class T>
  static constexpr T INITIAL_CONSTANT_WITH_PRECISION();

  template <>
  constexpr double INITIAL_CONSTANT_WITH_PRECISION<double>() {return std::numeric_limits<double>::max() / 16;}

  template <>
  constexpr float INITIAL_CONSTANT_WITH_PRECISION<float>() {return std::numeric_limits<float>::max() / 16;}

  template <class T>
  static constexpr T MIN_ACCEPTED_WITH_PRECISION();

  template <>
  constexpr double MIN_ACCEPTED_WITH_PRECISION<double>() {return 0.l;}

  template <>
  constexpr float MIN_ACCEPTED_WITH_PRECISION<float>() {return 1e-28f;}
}

template <class PRECISION>
class PairhmmImpl {

 public:
  size_t max_original_read_length = 0;  // updated during calculate()
  size_t max_padded_read_length = 0;    // updated during calculate()

  PairhmmImpl(const size_t initial_size = INITIAL_SIZE) :
    constants {initial_size},
    diagonals {initial_size},
    ph2pr{} {
    for (auto i=static_cast<PRECISION>(0); i!=MAX_PH2PR_INDEX; ++i) 
      ph2pr.push_back(pow(static_cast<PRECISION>(10.0), (-i) / static_cast<PRECISION>(10.0))); 
  }

  std::vector<double> calculate (const Testcase& testcase) {
    max_original_read_length = calculate_max_read_length(testcase.reads);
    const auto padded_haplotypes = pad_haplotypes(testcase.haplotypes);
    const auto padded_reads = pad_reads(testcase.reads);
    max_padded_read_length = calculate_max_read_length(padded_reads);
    resize_constants(max_padded_read_length);
    diagonals.resize(max_padded_read_length);
    return calculate(padded_reads, padded_haplotypes);
  }

  void recalculate (const Testcase& testcase, std::vector<double>& results, const size_t max_original_read_length_, const size_t max_padded_read_length_) {
    max_original_read_length = max_original_read_length_;
    max_padded_read_length = max_padded_read_length_;
    auto padded_pairs = std::vector<IndexedPair>{};
    padded_pairs.reserve(results.size());
    auto padded_read_cache = std::unordered_map<uint32_t, Read<PRECISION>>{};
    padded_read_cache.reserve(testcase.reads.size());
    auto padded_haplotype_cache = std::unordered_map<uint32_t, Haplotype>{};
    padded_haplotype_cache.reserve(testcase.haplotypes.size());
    auto master_idx = 0u;
    for (auto r : results) {
      if (r == FAILED_RUN_RESULT) {
        const auto read_idx = master_idx / testcase.haplotypes.size();
        const auto hap_idx = master_idx % testcase.haplotypes.size();
        const auto& padded_read = from_cache(padded_read_cache, read_idx, testcase.reads[read_idx]);
        const auto& padded_haplotype = from_cache(padded_haplotype_cache, hap_idx, testcase.haplotypes[hap_idx]);
        padded_pairs.emplace_back(master_idx, &padded_read, &padded_haplotype);
      }
      ++master_idx;
    }
    recalculate(padded_pairs, results, testcase.haplotypes.size());
  }

 protected:
  std::vector<Constants<PRECISION>> constants;
  Diagonals<PRECISION> diagonals;
  std::vector<PRECISION> ph2pr;

  static constexpr auto MAX_PH2PR_INDEX = 128;
  static constexpr auto INITIAL_CONSTANT = constants_with_precision::INITIAL_CONSTANT_WITH_PRECISION<PRECISION>();
  static constexpr auto INITIAL_SIZE = 250;
  static constexpr auto MIN_ACCEPTED = constants_with_precision::MIN_ACCEPTED_WITH_PRECISION<PRECISION>();  
  static constexpr auto FAILED_RUN_RESULT = std::numeric_limits<double>::min();

  struct IndexedPair {
    const uint32_t index;
    const Read<PRECISION> * const read;
    const Haplotype * const haplotype;

    IndexedPair(const uint32_t index_, const Read<PRECISION> * const read_, const Haplotype * const haplotype_) :
      index{index_},
      read{read_},
      haplotype{haplotype_}
    {}
  };

  template<class T>
  void pad (T& v, size_t padding) const {
    v.insert(v.end(), padding, 0);
  }

  Haplotype pad_and_reverse_haplotype(const Haplotype& haplotype, const size_t left_padding, const size_t right_padding) const {
    const auto padded_length = left_padding + haplotype.bases.size() + right_padding;
    auto p = Haplotype{};
    p.original_length = haplotype.bases.size();  
    p.bases.reserve(padded_length);
    pad(p.bases, left_padding);
    p.bases.insert(p.bases.end(), haplotype.bases.rbegin(), haplotype.bases.rend());
    pad(p.bases, right_padding);
    return p;
  }

  const Haplotype pad_haplotype(const Haplotype& haplotype) const {
    const auto left_padding = max_original_read_length + 1;
    const auto right_padding = max_original_read_length;
    return pad_and_reverse_haplotype(haplotype, left_padding, right_padding);
  }

  template<class Result_Type, const bool convert = true>
  std::vector<Result_Type> pad_and_convert_qual(const std::vector<uint8_t>& qual, const size_t left_padding, const size_t right_padding = 0) const {
    const auto padded_length = left_padding + qual.size() + right_padding;
    auto p = std::vector<Result_Type>{};
    p.reserve(padded_length);
    pad(p, left_padding);
    for (const auto q : qual)
      p.push_back(convert ? ph2pr[q] : q);
    pad(p, right_padding);
    return p;
  }

  Read<PRECISION> pad_read(const Read<uint8_t>& read) const {
    const auto left_padding = 1;
    auto padded_read = Read<PRECISION>{};
    padded_read.original_length = read.bases.size();
    padded_read.bases      = pad_and_convert_qual<uint8_t, false>(read.bases, left_padding);
    padded_read.base_quals = pad_and_convert_qual<PRECISION>(read.base_quals, left_padding);
    padded_read.ins_quals  = pad_and_convert_qual<PRECISION>(read.ins_quals, left_padding);
    padded_read.del_quals  = pad_and_convert_qual<PRECISION>(read.del_quals, left_padding);
    padded_read.gcp_quals  = pad_and_convert_qual<PRECISION>(read.gcp_quals, left_padding);
    return padded_read;
  }

  std::vector<Read<PRECISION>> pad_reads(const std::vector<Read<uint8_t>>& reads) const {
    auto padded_reads = std::vector<Read<PRECISION>>{};
    padded_reads.reserve(reads.size());
    for (auto& read : reads) 
      padded_reads.push_back(pad_read(read));
    return padded_reads;
  }

  std::vector<Haplotype> pad_haplotypes(const std::vector<Haplotype>& haplotypes) const {
    auto padded_haplotypes = std::vector<Haplotype>{};
    padded_haplotypes.reserve(haplotypes.size());
    for (const auto& haplotype : haplotypes) 
      padded_haplotypes.push_back(pad_haplotype(haplotype));
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

  void update_constants(const Read<PRECISION>& read) {
    const auto rows = read.bases.size();
    for (auto i = 1u; i != rows; ++i) {
      constants[i].mm = static_cast<PRECISION>(1) - (read.ins_quals[i] + read.del_quals[i]);
      constants[i].gm = static_cast<PRECISION>(1) - read.gcp_quals[i];
      constants[i].mx = read.ins_quals[i];
      constants[i].xx = read.gcp_quals[i];
      constants[i].my = i < (rows - 1) ? read.del_quals[i] : static_cast<PRECISION>(1);
      constants[i].yy = i < (rows - 1) ? read.gcp_quals[i] : static_cast<PRECISION>(1);
    }
  }

  void update_diagonals(const Haplotype& hap) {
    diagonals.zero();
    const auto first_row_value = INITIAL_CONSTANT / hap.original_length;
    diagonals.ypp[0] = diagonals.yp[0] = diagonals.y[0] = first_row_value;
  }

  inline PRECISION calculate_prior(const char read_base, const char hap_base, const PRECISION base_qual) const {
    return  ((read_base == hap_base) || (read_base == 'N') || (hap_base == 'N')) ? 
      static_cast<PRECISION>(1) - base_qual : 
      base_qual;
  }

  double compute_full_prob(const Read<PRECISION>& read, const Haplotype& haplotype) {
    const auto rows = read.bases.size();        // number of rows in the diagonals (padded read length)
    const auto hl = haplotype.original_length;  // haplotype original length (unpadded)
    const auto rl = read.original_length;       // read original length (unpadded)
    const auto mrl = max_original_read_length;  // alias for max original read length for readability in the code below (max read length in the testcase)
    const auto fd = mrl - rl;                   // first diagonal to compute (saves compute of all-0 diagonals when read is shorter than the padding - which will be the maximum read length in the testcase)
    auto result = 0.l;                          // result accumulator
    for (auto d = fd; d != mrl + hl - 1; ++d) { // d for diagonal
      for (auto r = 1u; r != rows; ++r) {       // r for row
        const auto hap_idx = mrl+hl-d+r-1;
        const auto read_base = read.bases[r];
        const auto hap_base = haplotype.bases[hap_idx];
        const auto prior = calculate_prior(read_base, hap_base, read.base_quals[r]);
        diagonals.m[r] = prior * ((diagonals.mpp[r-1] * constants[r].mm) + (constants[r].gm * (diagonals.xpp[r-1] + diagonals.ypp[r-1])));
        diagonals.x[r] = diagonals.mp[r-1] * constants[r].mx + diagonals.xp[r-1] * constants[r].xx; 
        diagonals.y[r] = diagonals.mp[r] * constants[r].my + diagonals.yp[r] * constants[r].yy; 
      }
      result += diagonals.m[rows-1] + diagonals.x[rows-1];
      diagonals.rotate();
    }
    return result < MIN_ACCEPTED ? 
      FAILED_RUN_RESULT :                       // if we underflowed return failed constant to rerun with higher precision if desired
      log10(static_cast<double>(result)) - log10(static_cast<double>(INITIAL_CONSTANT));
  }

  std::vector<double> calculate (const std::vector<Read<PRECISION>>& padded_reads, const std::vector<Haplotype>& padded_haplotypes) {
    auto results = std::vector<double>{};
    results.reserve(padded_reads.size() * padded_haplotypes.size());
    for (const auto& read : padded_reads) {
      update_constants(read);
      for (const auto& hap : padded_haplotypes) {
        update_diagonals(hap);
        results.push_back(compute_full_prob(read, hap));
      }
    }
    return results;
  }

  void recalculate (const std::vector<IndexedPair>& padded_pairs, std::vector<double>& results, const size_t n_haplotypes) {
    auto last_read_idx = -1;
    for (const auto& pp : padded_pairs) {
      auto read_idx = pp.index / n_haplotypes;
      if (read_idx != last_read_idx) 
        update_constants(*pp.read);    // only update if this is a new read (reads come in sorted)! 
      update_diagonals(*pp.haplotype); // should always update
      results[pp.index] = compute_full_prob(*pp.read, *pp.haplotype);
      last_read_idx = read_idx;
    }
  }

  const Read<PRECISION>& from_cache(std::unordered_map<uint32_t, Read<PRECISION>>& cache, const uint32_t index, const Read<uint8_t>& original) {
    if (cache.count(index) == 0) 
      cache.emplace(index, pad_read(original));
    return cache.at(index);
  }

  const Haplotype& from_cache(std::unordered_map<uint32_t, Haplotype>& cache, const uint32_t index, const Haplotype& original) {
    if (cache.count(index) == 0) 
      cache.emplace(index, pad_haplotype(original));
    return cache.at(index);
  }

};

#endif

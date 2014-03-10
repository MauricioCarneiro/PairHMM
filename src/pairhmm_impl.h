#ifndef PAIRHMM_IMPL_H
#define PAIRHMM_IMPL_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <ostream>
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

template <class PRECISION, class DIAGONALS, class CONSTANTS, int VECSIZE = 1>
class PairhmmImpl {
  size_t m_max_original_read_length = 0;  // updated during calculate()
  size_t m_max_padded_read_length = 0;    // updated during calculate()

public:
  size_t max_padded_read_length(void) const { return m_max_padded_read_length; }
  size_t max_original_read_length(void) const { return m_max_original_read_length; }

  void set_max_padded_read_length(size_t length) const { m_max_padded_read_length = length; }
  void set_max_original_read_length(size_t length) const { m_max_original_read_length = length; }

  const std::vector<Haplotype<PRECISION>>& padded_haplotypes() const { return m_padded_haplotypes; }

  PairhmmImpl(const size_t initial_size = INITIAL_SIZE) :
    m_constants {initial_size},
    m_diagonals {initial_size},
    m_ph2pr{} {
    m_ph2pr.reserve(MAX_PH2PR_INDEX);
    for (auto i=static_cast<PRECISION>(0); i!=MAX_PH2PR_INDEX; ++i)
      m_ph2pr.push_back(pow(static_cast<PRECISION>(10.0), (-i) / static_cast<PRECISION>(10.0)));
  }

  virtual ~PairhmmImpl() { }

  std::vector<double> calculate (const Testcase& testcase) {
    m_max_original_read_length = calculate_max_read_length(testcase.reads);
    m_padded_haplotypes = pad_haplotypes(testcase.haplotypes);
    const auto padded_reads = pad_reads(testcase.reads);
    m_max_padded_read_length = calculate_max_read_length(padded_reads);
    m_constants.reserve(m_max_padded_read_length);
    m_diagonals.reserve(m_max_padded_read_length);
    return calculate(padded_reads);
  }

  void recalculate (const Testcase& testcase, std::vector<double>& results, const size_t max_original_read_length, 
      const size_t max_padded_read_length) {
    m_max_original_read_length = max_original_read_length;
    m_max_padded_read_length = max_padded_read_length;
    auto master_idx = 0;
    auto padded_read = Read<PRECISION,PRECISION>{};
    for (auto read : testcase.reads) {
      auto has_padded_read = false;
      for (auto hap_idx = 0u; hap_idx != testcase.haplotypes.size(); ++hap_idx) {
        if (results[master_idx] == FAILED_RUN_RESULT) {
          if (!has_padded_read) {
            padded_read = pad_read(read); // making a copy here so I can keep it across haplotype iterations
            m_constants.update(padded_read);
            has_padded_read = true;
          }
          auto padded_haplotype = pad_haplotype(testcase.haplotypes[hap_idx]); // can try to optimize this later, but only after benchmarking
          m_diagonals.update(INITIAL_CONSTANT/padded_haplotype.original_length);
          results[master_idx] = compute_full_prob(padded_read, padded_haplotype);
        }
        ++master_idx;
      }
    }
  }

 protected:
  CONSTANTS m_constants;
  DIAGONALS m_diagonals;
  std::vector<PRECISION> m_ph2pr;
  std::vector<Haplotype<PRECISION>> m_padded_haplotypes;

  static constexpr auto MAX_PH2PR_INDEX = 128;
  static constexpr auto INITIAL_CONSTANT = constants_with_precision::INITIAL_CONSTANT_WITH_PRECISION<PRECISION>();
  static constexpr auto INITIAL_SIZE = 250;
  static constexpr auto MIN_ACCEPTED = constants_with_precision::MIN_ACCEPTED_WITH_PRECISION<PRECISION>();
  static constexpr auto FAILED_RUN_RESULT = std::numeric_limits<double>::min();

  template<class T>
    void pad (T& v, size_t padding) const {
      v.insert(v.end(), padding, 0);
    }

  Haplotype<PRECISION> pad_and_reverse_haplotype(const Haplotype<uint8_t>& haplotype, const size_t left_padding, const size_t right_padding) const {
    const auto padded_length = left_padding + haplotype.bases.size() + right_padding;
    auto p = Haplotype<PRECISION>{};
    p.original_length = haplotype.bases.size();
    p.bases.reserve(padded_length);
    pad(p.bases, left_padding);
    p.bases.insert(p.bases.end(), haplotype.bases.rbegin(), haplotype.bases.rend());
    pad(p.bases, right_padding);
    return p;
  }

  const Haplotype<PRECISION> pad_haplotype(const Haplotype<uint8_t>& haplotype) const {
    const auto left_padding = m_max_original_read_length + 1;
    const auto right_padding = m_max_original_read_length + (VECSIZE-1);
    return pad_and_reverse_haplotype(haplotype, left_padding, right_padding);
  }

  template<class Result_Type, const bool convert = true>
    std::vector<Result_Type> pad_and_convert_qual(const std::vector<uint8_t>& qual, const size_t left_padding, const size_t right_padding) const {
      const auto padded_length = left_padding + qual.size() + right_padding;
      auto p = std::vector<Result_Type>{};
      p.reserve(padded_length);
      pad(p, left_padding);
      for (const auto q : qual)
        p.push_back(convert ? m_ph2pr[q] : q);
      pad(p, right_padding);
      return p;
    }

  Read<PRECISION,PRECISION> pad_read(const Read<uint8_t,uint8_t>& read) const {
    const auto left_padding = 1;
    const auto right_padding = VECSIZE-1; // enough for worst-case
    auto padded_read = Read<PRECISION,PRECISION>{};
    padded_read.left_padding = left_padding;
    padded_read.original_length = read.bases.size();
    padded_read.right_padding = right_padding;
    padded_read.bases      = pad_and_convert_qual<PRECISION, false>(read.bases, left_padding, right_padding);
    padded_read.base_quals = pad_and_convert_qual<PRECISION>(read.base_quals, left_padding, right_padding);
    padded_read.ins_quals  = pad_and_convert_qual<PRECISION>(read.ins_quals, left_padding, right_padding);
    padded_read.del_quals  = pad_and_convert_qual<PRECISION>(read.del_quals, left_padding, right_padding);
    padded_read.gcp_quals  = pad_and_convert_qual<PRECISION>(read.gcp_quals, left_padding, right_padding);
    return padded_read;
  }

  std::vector<Read<PRECISION,PRECISION>> pad_reads(const std::vector<Read<uint8_t, uint8_t>>& reads) const {
    auto padded_reads = std::vector<Read<PRECISION,PRECISION>>{};
    padded_reads.reserve(reads.size());
    for (auto& read : reads)
      padded_reads.push_back(pad_read(read));
    return padded_reads;
  }

  std::vector<Haplotype<PRECISION>> pad_haplotypes(const std::vector<Haplotype<uint8_t>>& haplotypes) const {
    auto padded_haplotypes = std::vector<Haplotype<PRECISION>>{};
    padded_haplotypes.reserve(haplotypes.size());
    for (const auto& haplotype : haplotypes)
      padded_haplotypes.push_back(pad_haplotype(haplotype));
    return padded_haplotypes;
  }

  template <class T, class U>
    size_t calculate_max_read_length(const std::vector<Read<T,U>>& reads) const {
      auto max_read_length = 0u;
      for (const auto& read : reads)
        max_read_length = max_read_length >= read.bases.size() ? max_read_length : read.bases.size();
      return max_read_length;
    }

  double compute_full_prob(const Read<PRECISION,PRECISION>& read, const Haplotype<PRECISION>& haplotype) {
    return do_compute_full_prob(read, haplotype);
  }

  std::vector<double> calculate (const std::vector<Read<PRECISION,PRECISION>>& padded_reads) {
    auto results = std::vector<double>{};
    results.reserve(padded_reads.size() * m_padded_haplotypes.size());
    for (const auto& read : padded_reads) {
      m_constants.update(read);
      for (const auto& hap : m_padded_haplotypes) {
        m_diagonals.update(INITIAL_CONSTANT / hap.original_length);
        results.push_back(compute_full_prob(read, hap));
      }
    }
    return results;
  }

protected:

  virtual double do_compute_full_prob(const Read<PRECISION,PRECISION>& read, const Haplotype<PRECISION>& haplotype) = 0;

};

#endif

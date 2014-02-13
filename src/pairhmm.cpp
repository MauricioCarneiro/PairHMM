#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

#include "pairhmm.h"
#include "testcase.h"
#include "read.h"
#include "haplotype.h"
#include "constants.h"
#include "diagonals.h"

using namespace std;

template<typename T>
vector<T> pad(const vector<T>& v, const uint32_t leading_pad = 1, const uint32_t trailing_pad = 0) {
  auto r = vector<T> {};
  r.resize(v.size() + leading_pad + trailing_pad);
  copy(v.begin(), v.end(), r.begin() + leading_pad);
  return r;
}

vector<float> real(const vector<uint8_t>& v, const vector<float>& ph2pr) {
  auto result = vector<float>{};
  for (auto& q :v)
    result.push_back(ph2pr[q]);
  return result;
}

vector<uint8_t> invert(const vector<uint8_t>& v) {
  auto r = vector<uint8_t>{};
  for(auto element = v.rbegin(); element != v.rend(); ++element)
    r.push_back(*element);
  return r;
}

Haplotype pad(const Haplotype& haplotype, const uint32_t read_length) {
  return Haplotype{pad(invert(haplotype.get_bases()), read_length + 1, read_length), haplotype.get_padded_length()};
}

Read<float> pad(const Read<uint8_t>& read, const vector<float>& php2pr) {
  auto b = pad(read.get_bases());
  auto q = pad(real(read.get_base_quals(), php2pr));
  auto i = pad(real(read.get_ins_quals(), php2pr));
  auto d = pad(real(read.get_del_quals(), php2pr));
  auto g = pad(real(read.get_gcp_quals(), php2pr));
  return Read<float>{b,q,i,d,g,read.get_padded_length()};
}

vector<Haplotype> pad_haplotypes(const vector<Haplotype>& haplotypes, const size_t max_read_length) {
  auto padded_haplotypes = vector<Haplotype>{};
  for (const auto& haplotype : haplotypes) 
    padded_haplotypes.push_back(pad(haplotype, max_read_length));
  return padded_haplotypes;
}

template<class precision>
Pairhmm<precision>::Pairhmm(const size_t initial_size) :
    constants {initial_size},
    diagonals {initial_size},
    ph2pr{}
  {
    for (auto i=0.f; i!=MAX_VAL; ++i) 
      ph2pr.push_back(pow(static_cast<precision>(10.0), (-i) / static_cast<precision>(10.0))); 
  }

template<class precision>
double Pairhmm<precision>::compute_full_prob(const Read<float>& read, const Haplotype& haplotype) {
  const auto rows = read.get_padded_length();
  const auto rl = read.get_original_length();
  const auto hl = haplotype.get_original_length();
  const auto first_row_value = INITIAL_CONSTANT / haplotype.get_original_length();

  // initialize constant arrays
  if (constants.capacity() < rows) 
    constants.resize(rows);

  for (auto i = 1u; i != rows; ++i) {
    constants[i].mm = 1.f - (read.get_ins_quals()[i] + read.get_del_quals()[i]);
    constants[i].gm = 1.f - read.get_gcp_quals()[i];
    constants[i].mx = read.get_ins_quals()[i];
    constants[i].xx = read.get_gcp_quals()[i];
    constants[i].my = i < (rows - 1) ? read.get_del_quals()[i] : 1.f;
    constants[i].yy = i < (rows - 1) ? read.get_gcp_quals()[i] : 1.f;
  }

  // initialize hmm matrix (diagonal vectors)
  auto m   = vector<float>(rows,0);
  auto mp  = vector<float>(rows,0);
  auto mpp = vector<float>(rows,0);
  auto x   = vector<float>(rows,0);
  auto xp  = vector<float>(rows,0);
  auto xpp = vector<float>(rows,0);
  auto y   = vector<float>(rows,0);
  auto yp  = vector<float>(rows,0);
  auto ypp = vector<float>(rows,0);
  ypp[0] = yp[0] = y[0] = first_row_value;

  // main loop
  auto result = 0.l;
  for (auto d = 0u; d != rl + hl - 1; ++d) {  // d for diagonal
    m[0] = 0.f;
    x[0] = 0.f;
    y[0] = first_row_value;
    for (auto r = 1u; r != rows; ++r) {  // r for row
      auto hap_idx = rl+hl-d+r-1;
      auto read_base = read.get_bases()[r];
      auto hap_base = haplotype.get_bases()[hap_idx];
      auto distm = ( (read_base == hap_base) || (read_base == 'N') || (hap_base == 'N') ) ? 1.f - read.get_base_quals()[r] : read.get_base_quals()[r];
      m[r] = distm * ((mpp[r-1] * constants[r].mm) + (constants[r].gm * (xpp[r-1] + ypp[r-1])));
      x[r] = mp[r-1] * constants[r].mx + xp[r-1] * constants[r].xx; 
      y[r] = mp[r] * constants[r].my + yp[r] * constants[r].yy; 
    }
    result += m[rows-1] + x[rows-1];
    swap(mpp, mp);
    swap(xpp, xp);
    swap(ypp, yp);
    swap(mp, m); 
    swap(xp, x);
    swap(yp, y);
  }
  return log10(static_cast<double>(result)) - log10(static_cast<double>(INITIAL_CONSTANT));
}

template<class precision>
vector<Read<float>> Pairhmm<precision>::pad_reads(const vector<Read<uint8_t>>& reads) {
  auto padded_reads = vector<Read<float>>{};
  for (const auto& read : reads) 
    padded_reads.push_back(pad(read, ph2pr));
  return padded_reads;
}

template<class precision>
vector<double> Pairhmm<precision>::calculate (const vector<Read<float>>& padded_reads, const vector<Haplotype>& padded_haplotypes) {
  auto results = vector<double>{};
  for (const auto& hap : padded_haplotypes) 
    for (const auto& read : padded_reads) 
      results.push_back(compute_full_prob(read, hap, constants));
  return results;
}

template<class precision>
vector<double> Pairhmm<precision>::calculate (Testcase& testcase) {
  const auto max_read_length = testcase.calculate_max_read_length();
  const auto padded_haplotypes = pad_haplotypes(testcase.get_haplotypes(), max_read_length);
  const auto padded_reads = pad_reads(testcase.get_reads());
  return pairhmm(padded_reads, padded_haplotypes);
}


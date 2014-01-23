#include <istream>
#include <string>
#include <vector>
#include <limits>

#include "testcase_iterator.h"
#include "read.h"
#include "haplotype.h"

using namespace std;

TestcaseIterator::TestcaseIterator(istream* const input) :
  input_stream { input },
  current      { fetch_next_element() }
{}

TestcaseIterator::TestcaseIterator(TestcaseIterator&& original) :
  input_stream { original.input_stream },
  current {move(original.current) }
{
  original.input_stream = nullptr;
}

bool TestcaseIterator::operator!=(const TestcaseIterator& rhs) {
  return input_stream != rhs.input_stream;
}

Testcase& TestcaseIterator::operator*() {
  return current;
}

Testcase& TestcaseIterator::operator++() {
  current = fetch_next_element();
  return current;
}

Testcase TestcaseIterator::fetch_next_element() {
  const auto test_size = parse_test_size();
  auto reads = vector<Read>{};
  auto haplotypes = vector<Haplotype>{};

  for (size_t i = 0; i < test_size.first; ++i)  // parse all reads
    reads.push_back(parse_read()); 
  for (size_t i = 0; i < test_size.second; ++i) // parse all haplotypes
    haplotypes.push_back(parse_haplotype());

  return Testcase {move(reads), move(haplotypes)};
}

pair<size_t, size_t> TestcaseIterator::parse_test_size() {
  int n_reads, n_haplotypes;
  *input_stream >> n_reads; 
  *input_stream >> n_haplotypes; 
  return make_pair(n_reads, n_haplotypes);
}

Read TestcaseIterator::parse_read() {
  const auto read_bases      = parse_string();
  const auto read_quals      = parse_string();
	const auto read_ins_quals  = parse_string();
	const auto read_del_quals  = parse_string();
	const auto read_gcp_quals  = parse_string();
  return Read{read_bases, read_quals, read_ins_quals, read_del_quals, read_gcp_quals};
}

Haplotype TestcaseIterator::parse_haplotype() {
	const auto haplotype_bases = parse_string();
  return Haplotype {haplotype_bases};
}

string TestcaseIterator::parse_string() {
  string s;
  *input_stream >> s;
  return s;
}



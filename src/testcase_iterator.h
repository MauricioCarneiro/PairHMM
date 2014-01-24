#ifndef __TESTCASE_ITERATOR__
#define __TESTCASE_ITERATOR__

#include <istream>
#include "testcase.h"

class TestcaseIterator {
  std::istream * input_stream; ///< a pointer to the input stream
  Testcase current;

  Testcase fetch_next_element();
  std::pair<size_t, size_t> parse_test_size();
  Read parse_read();
  Haplotype parse_haplotype();
  std::string parse_string();

 public:
  
  TestcaseIterator() :
    input_stream {nullptr},
    current {}
  {}

  TestcaseIterator(std::istream* const input);
  TestcaseIterator(const TestcaseIterator&) = delete;
  TestcaseIterator(TestcaseIterator&& original);

  bool operator!=(const TestcaseIterator&);
  Testcase& operator*();
  Testcase& operator++();

};

#endif

#ifndef __SA_H__
#define __SA_H__

#include <iostream>
#include <string>
#include <vector>

class SuffixArray {
 private:
  const std::string ref;
  std::vector<size_t> suffixes;
  std::vector<size_t> lcp;
  int isPrefix(const std::string &mot, size_t pos) const;
  size_t longest_common_prefix(size_t p1, size_t p2) const;
  void init();
 public:
  SuffixArray(const std::string &txt);
  std::string at(size_t i) const;
  std::string operator[](size_t i) const;
  size_t lookup(const std::string &mot) const;
  size_t nb_occurrences(const std::string &mot) const;
  void toStream(std::ostream &os) const;
  size_t size() const;
};

std::ostream &operator<<(std::ostream &os, const SuffixArray &sa);

#endif

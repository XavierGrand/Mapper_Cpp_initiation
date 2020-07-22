#include "SA_ES.h"
#include <algorithm>
#include <libgen.h>

using namespace std;

struct sa_compare_pos {
  const EncodedSequence &ref;
  sa_compare_pos(const EncodedSequence &ref):ref(ref) {
  }
  bool operator()(size_t a, size_t b) {
    int v = 0;
    while (!v && (a < ref.size()) && (b < ref.size())) {
      if (ref[a] < ref[b]) {
        v = -1;
      } else {
        if (ref[a] > ref[b]) {
          v = 1;
        } else {
          ++a;
          ++b;
        }
      }
    }
    if (!v) {
      v = (a > b) ? -1 : 1;
    }
    return (v < 0);
  }
};

int SuffixArray::isPrefix(const string &mot, size_t pos) const {
  size_t p = 0;
  int res = 0;
  while (!res && (p < mot.length())) {
    if (mot[p] < ref[pos + p]) {
      res = -1;
    } else {
      if (mot[p] > ref[pos + p]) {
        res = 1;
      }
    }
    ++p;
  }
  return res;
}

size_t SuffixArray::lookup(const string &mot) const {
  size_t deb = 0, fin = suffixes.size();
  bool found = false;
  while (!found && (deb < fin)) {
    size_t p = (deb + fin) / 2;
    int res = isPrefix(mot, suffixes[p]);
    if (!res) {
      deb = fin = p;
      found = true;
    } else {
      if (res == 1) {
        deb = p + 1;
      } else {
        fin = p;
      }
    }
  }
  return (found ? deb : (size_t) -1);
}

size_t SuffixArray::nb_occurrences(const string &mot) const {
  size_t p = lookup(mot);
  if (p == (size_t) -1) {
    return 0;
  }
  size_t cpt = 1, i = p;
  
  // Remonter à partir de p vers le haut tant que mot est préfixe du suffixe començant en SA[p]
  while (i && (lcp[i] >= mot.length())) {
    ++cpt;
    --i;
  }
  // Descendre à partir de p vers le bas tant que mot est préfixe du suffixe començant en SA[p]
  i = p + 1;
  while ((i < lcp.size()) && (lcp[i] >= mot.length())) {
    ++cpt;
    ++i;
  }
  return cpt;
}

void SuffixArray::toStream(ostream &os) const {
  os << "SA :\n";
  for (size_t i = 0; i < suffixes.size(); ++i) {
    os << "- " << i << ": [" << suffixes[i] << "] (='";
    for (size_t j = suffixes[i]; j < ref.size(); ++j) {
      os << ref[j];
    }
    os << "')"
       << " [LCP = " << lcp[i] << "]\n";
  }
  os << endl;
}

size_t SuffixArray::longest_common_prefix(size_t p1, size_t p2) const {
  size_t cpt = 0;
  bool ok = true;
  while (ok && ((p1 < ref.size()) && (p2 < ref.size()))) {
    ok = (ref[p1] == ref[p2]);
    cpt += ok;
    ++p1;
    ++p2;
  }
  return cpt;
}

void SuffixArray::init() {
  suffixes.resize(ref.size());
  lcp.resize(ref.size());
  for (size_t i = 0; i < suffixes.size(); ++i) {
    suffixes[i] = i;
  }
  sort(suffixes.begin(), suffixes.end(), sa_compare_pos(ref));
  lcp[0] = 0;
  for (size_t i = 1; i < lcp.size(); ++i) {
    lcp[i] = longest_common_prefix(suffixes[i - 1], suffixes[i]);
  }  
}

SuffixArray::SuffixArray(const EncodedSequence &encseq):ref(encseq) {
  init();
}

ostream &operator<<(ostream &os, const SuffixArray &sa) {
  sa.toStream(os);
  return os;
}


#ifdef TEST_SA_ES_MAIN
int main(int argc, char **argv) {
  if (argc != 2) {
    cerr << "usage: " << basename(argv[0]) << " <DNA_Seq>" << endl;
    return 1;
  }
  
  SuffixArray SA(argv[1]);
  cout << SA << endl;
  cout << "Recherche de ACA => " << SA.lookup("ACA") << endl;
  cout << "Nb occurrences de ACA => " << SA.nb_occurrences("ACA") << endl;
  cout << "Recherche de CT => " << SA.lookup("CT") << endl;
  cout << "Nb occurrences de CT => " << SA.nb_occurrences("CT") << endl;
  cout << "Recherche de GCT => " << SA.lookup("GCT") << endl;
  cout << "Nb occurrences de GCT => " << SA.nb_occurrences("GCT") << endl;
  return 0;
}

#endif

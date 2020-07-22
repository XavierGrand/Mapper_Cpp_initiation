#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

void PrintSA_naive(const vector<string> &sa) {
  cout << "SA :\n";
  for (size_t i = 0; i < sa.size(); ++i) {
    cout << "- " << i << ": '" << sa[i] << "'\n";
  }
  cout << endl;
}

int isPrefix_naif(const string &mot1, const string &mot2) {
  size_t p = 0;
  int res = 0;
  while (!res && (p < mot1.length())) {
    if (mot1[p] < mot2[p]) {
      res = -1;
    } else {
      if (mot1[p] > mot2[p]) {
        res = 1;
      }
    }
    ++p;
  }
  return res;
}

bool lookup_naif(const vector<string> &sa, const string &mot) {
  size_t deb = 0, fin = sa.size();
  bool found = false;
  while (!found && (deb < fin)) {
    size_t p = (deb + fin) / 2;
    int res = isPrefix_naif(mot, sa[p]);
    if (!res) {
      found = true;
    } else {
      if (res == 1) {
        deb = p + 1;
      } else {
        fin = p;
      }
    }
  }
  return found;
}


vector<string> SA_naive(const string &txt) {
  
  vector<string> SA(txt.length());

  for (size_t i = 0; i < SA.size(); ++i) {
    SA[i] = txt.substr(i);
  }

  PrintSA_naive(SA);

  sort(SA.begin(), SA.end());

  cout << "Après tri:" << endl;

  PrintSA_naive(SA);
  return SA;
}


struct sa_compare_pos {
  const string &ref;
  sa_compare_pos(const string &ref):ref(ref) {
    cerr << "Création d'un objet de comparaison pour le texte '" << ref << "'" << endl;
  }
  bool operator()(size_t a, size_t b) {
    cerr << "L'opérateur de comparaison du texte '" << ref
         << "' est appelé avec a = " << a
         << " et b = " << b << endl;
    int v = 0;
    while (!v && (a < ref.length()) && (b < ref.length())) {
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

int isPrefix_subtil(const string &mot1, const string &mot2, size_t pos) {
  size_t p = 0;
  int res = 0;
  while (!res && (p < mot1.length())) {
    if (mot1[p] < mot2[pos + p]) {
      res = -1;
    } else {
      if (mot1[p] > mot2[pos + p]) {
        res = 1;
      }
    }
    ++p;
  }
  return res;
}

bool lookup_subtil(const vector<size_t> &sa, const string &ref, const string &mot) {
  size_t deb = 0, fin = sa.size();
  bool found = false;
  while (!found && (deb < fin)) {
    size_t p = (deb + fin) / 2;
    int res = isPrefix_subtil(mot, ref, sa[p]);
    if (!res) {
      found = true;
    } else {
      if (res == 1) {
        deb = p + 1;
      } else {
        fin = p;
      }
    }
  }
  return found;
}
void PrintSA_subtile(const vector<size_t> &sa, const string &ref) {
  cout << "SA :\n";
  for (size_t i = 0; i < sa.size(); ++i) {
    cout << "- " << i << ": [" << sa[i] << "] (='" << ref.substr(sa[i]) << "')\n";
  }
  cout << endl;
}

vector<size_t> SA_subtile(const string &txt) {
  
  vector<size_t> SA(txt.length());

  for (size_t i = 0; i < SA.size(); ++i) {
    SA[i] = i;
  }

  PrintSA_subtile(SA, txt);

  sa_compare_pos comparateur(txt);
  sort(SA.begin(), SA.end(), comparateur);

  cout << "Après tri:" << endl;

  PrintSA_subtile(SA, txt);
  return SA;
}

#ifdef TEST_SA_CONCEP_MAIN
int main(int argc, char **argv) {
  vector<string> SA_n = SA_naive("ROUDOUDOU");
  cout << "Recherche de ROU => " << lookup_naif(SA_n, "ROU") << endl;
  cout << "Recherche de OUX => " << lookup_naif(SA_n, "OUX") << endl;

  vector<size_t> SA_s = SA_subtile("ROUDOUDOU");
  cout << "Recherche de ROU => " << lookup_subtil(SA_s, "ROUDOUDOU", "ROU") << endl;
  cout << "Recherche de OUX => " << lookup_subtil(SA_s, "ROUDOUDOU", "OUX") << endl;
  return 0;
}

#endif

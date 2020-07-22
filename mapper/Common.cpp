#include "Common.h"
#include <string>

using namespace std;

bool Common::isNucl(char c) {
    return ((c == 'A') || (c == 'a') ||
            (c == 'C') || (c == 'c') ||
            (c == 'G') || (c == 'g') ||
            (c == 'T') || (c == 't') ||
            (c == 'U') || (c == 'u') ||
            (c == 'R') || (c == 'r') ||
            (c == 'Y') || (c == 'y') ||
            (c == 'S') || (c == 's') ||
            (c == 'W') || (c == 'w') ||
            (c == 'K') || (c == 'k') ||
            (c == 'M') || (c == 'm') ||
            (c == 'B') || (c == 'b') ||
            (c == 'D') || (c == 'd') ||
            (c == 'H') || (c == 'h') ||
            (c == 'V') || (c == 'v') ||
            (c == 'N') || (c == 'n') ||
            (c == '.') || (c == '-'));
}

bool Common::isAA(char c) {
    return ((c == 'A') || (c == 'a') ||
            (c == 'C') || (c == 'c') ||
            (c == 'D') || (c == 'd') ||
            (c == 'E') || (c == 'e') ||
            (c == 'F') || (c == 'f') ||
            (c == 'G') || (c == 'g') ||
            (c == 'H') || (c == 'h') ||
            (c == 'I') || (c == 'i') ||
            (c == 'K') || (c == 'k') ||
            (c == 'L') || (c == 'l') ||
            (c == 'M') || (c == 'm') ||
            (c == 'N') || (c == 'n') ||
            (c == 'P') || (c == 'p') ||
            (c == 'Q') || (c == 'q') ||
            (c == 'R') || (c == 'r') ||
            (c == 'S') || (c == 's') ||
            (c == 'T') || (c == 't') ||
            (c == 'V') || (c == 'v') ||
            (c == 'W') || (c == 'w') ||
            (c == 'Y') || (c == 'y') ||
            (c == 'B') || (c == 'b') ||
            (c == 'Z') || (c == 'z') ||
            (c == 'X') || (c == 'x') ||
            (c == '.') || (c == '-'));
}

string Common::reverseComplement(const string seq) { // const a ajouter, et il vaut mieux passer du flux que de copier la seq... passage par ref.
  size_t seqlen = seq.length();
  string revertedSeq = "";
  cout << "Sequence length = " << seqlen << endl;
  for (size_t i = seqlen; i > 0; i--) {
    if(seq[i] == 'A' || seq[i] == 'a') { //Que faire du cas des ARN ? En réalité, la retrotranscription d'un ARN en ARN n'a pas réellement d'intérêt, je ne crois pas que ce soit une réalité biologique (à vérifier, possible chez des virus), et s'il s'agit d'aligner sur une seq de ref, cette dernière sera de l'ADN... Donc T on laisse tomber le U.
      revertedSeq += "T";
    }
    else if(seq[i] == 'C' || seq[i] == 'c') {
      revertedSeq += "T";
    }
    else if(seq[i] == 'G' || seq[i] == 'g') {
      revertedSeq += "T";
    }
    else if(seq[i] == 'T' || seq[i] == 't' || seq[i] == 'U' || seq[i] == 'u') {
      revertedSeq += "A";
    }
    else{revertedSeq += "N";}
  }
  return revertedSeq;
}

char Common::revComp(const char* c) {
  char ch = *c;
  char res = 'N';
  if(ch == 'A' || ch == 'a') { //Que faire du cas des ARN ? En réalité, la retrotranscription d'un ARN en ARN n'a pas réellement d'intérêt, je ne crois pas que ce soit une réalité biologique (à vérifier, possible chez des virus), et s'il s'agit d'aligner sur une seq de ref, cette dernière sera de l'ADN... Donc T on laisse tomber le U.
    res = 'T';
  }
  else if(ch == 'C' || ch == 'c') {
    res = 'G';
  }
  else if(ch == 'G' || ch == 'g') {
    res = 'C';
  }
  else if(ch == 'T' || ch == 't' || ch == 'U' || ch == 'u') {
    res = 'A';
  }
  else{
    res = 'N';
  }
  return res;
}

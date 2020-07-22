#include "fastxfile.h"
#include "fastxseq.h"
#include "EncodedSequence.h"
#include "Common.h"
#include <vector>
#include <cassert>
#include <iostream>
#include <ios>

using namespace std;

void FastXFile::detectFormat() { //limpide
  format = UNKNOWN;
  bool detected = false;
  while (!detected && !ifs.eof()) {
    char s[256];
    ifs.read(s, 256);
    size_t nb = ifs.gcount();
    if (nb) {
      size_t i = 0;
      while (!detected && (i < nb)) {
        detected = (s[i] > ' ');//Dans la table ASCII, le caractère ' ' est de caractère qui précède les lettres, chiffre et autres... Avat, ce sont les caractère vides...
        ++i;
      }
      if (detected) {
        switch (s[--i]) {
        case ';':
        case '>': format = FASTA;
          break;
        case '@': format = FASTQ;
          break;
        default: cerr << "Warning: format non détecté, caractère '"
                      << s[i] << "' rencontré" << endl;
        }
      }
    }
  }
  if (format != UNKNOWN) {
    ifs.seekg(0);
  } else {
    ifs.close();
  }
}

const char* FastXFile::BioSeqFmt2String(BiologicalSequenceFormat fmt) {
  return (fmt == FASTA
          ? "FASTA"
          : (fmt == FASTQ
             ? "FASTQ"
             : "Unknown"));
}

FastXFile::FastXFile(const char *fname):fname(fname), ifs(), format(UNKNOWN) {
  if (fname) {
    open(fname);
  }
}

FastXFile::~FastXFile() {
  ifs.close();
}

void FastXFile::open(const char* fname) {
  if (fname) {
    ifs.close();
    ifs.open(fname);
    if (ifs) {
      detectFormat();
    } else {
      cerr << "Warning: Impossible d'ouvrir le fichier '" << fname << "'" << endl;
    }
  }
}

string FastXFile::getCurrentSequenceHeader() {
  streampos p = ifs.tellg();
  string entete;
  getline(ifs, entete);
  ifs.seekg(p);
  return entete;
}

void FastXFile::processCurrentSequence() { // FastXSeq // size_t, j'ai modifié la fonction pour qu'elle renvoie un objet FastXSeq.

  FastXSeq obj;

  //Les attributs pour construire l'objet FastXSeq :
  streampos entetePos(0);
  streampos seqPos(0);

  //Variables nécessaires au parsing :
  //attention si l'entete est trop long, plus de la taille du buffer, ça chie
  char buffer[4096]; 
  char previous_char = '\n';
  bool new_sequence = true;
  bool in_header = true;
  bool ok = true;
  
  while (ifs && ok) {
    ifs.read(buffer, 4096);
    size_t nb = ifs.gcount();
    size_t i = 0;
    if (nb) {
      // Virer les caractères vide en début de fichier. On pourra it mettre inférieur au '>' pour virer des éventuels commentaires qui auraient été ajoutés en amont des séquences.
      // while(buffer[i] < '!') { 
      //   ++i;
      // }
      //      clog << "Traitement du buffer = '" << buffer << "'" << endl;
      if (format == FASTA){
        if (new_sequence) {
          assert(buffer[0] == '>' || buffer[0] == ';');
          in_header = true;
          // Enregistrer l'adresse pour entete.
          entetePos = (streamoff)ifs.tellg() - (nb - i); 
        }
        while (in_header && (i < nb)) {
          previous_char = buffer[i];
          in_header = (buffer[i++] != '\n');
        }
        // clog << "Traitement du buffer à partir de la position " << i << " => '" << &buffer[i] << "'" << endl;
        if (!in_header) {
          new_sequence = false; 
          // On pourrait verifier que j'ai bien un caractère de l'alphabet IUPAC avant de recup l'adresse
          seqPos = (streamoff)ifs.tellg() - (nb - i);
          // la méthode renvoie true
          ok = !obj.parseA(ifs, entetePos, seqPos); // ok = !obj.parseA(ifs, entetePos, seqPos);
        }
      }
      else if (format == FASTQ){
        if (new_sequence) {
          // Début de l'entete si @ en début de ligne (problème avec @ c'est si une quality débute par @)
          assert(buffer[i] == '@' && previous_char == '\n'); 
          in_header = true;
          entetePos = (streamoff)ifs.tellg() - (nb - i);
        }
        while (in_header && (i < nb)) {
          previous_char = buffer[i];
          //On teste toute la ligne de header. On fait l'assomption que le header n'est que sur une ligne.
          in_header = (buffer[i++] != '\n'); 
        }
        if (!in_header) {
          new_sequence = false;
          seqPos = (streamoff)ifs.tellg() - (nb - i);
          ok = !obj.parseQ(ifs, entetePos, seqPos); // ok = !obj.parseQ(ifs, entetePos, seqPos);
        }
      }  
      else { cerr << "Unsupported format." << endl;
      }
    }
  }
  seqList.push_back(obj);
  // return obj;
}

ostream &operator<<(ostream &os, FastXFile::BiologicalSequenceFormat fmt) {
  return os << FastXFile::BioSeqFmt2String(fmt);
}

//Attention, il faut faire le lien avec la classe EncodedSequence pour lui laisser la responsabilité de la création d'un EncodedSequence.
vector<EncodedSequence> FastXFile::EncodedSequence(){
  vector<EncodedSequence> listES;
  // Passer toutes les FastXSeq du vecteur de FastXSeq créé à partir du fichier Fasta
  for (size_t u = 0; u < seqList.size(); u++) {
    // récupérer la position de la seq dans l'objet FastXSeq indice u du vecteur.
    char buffer[seqList[u].getNbNucl()];
    ifs.seekg(seqList[u].getSeqPos());
    ifs.read(buffer, seqList[u].getNbNucl());
    EncodedSequence es(buffer);
    listES.push_back(es);
  }
  return listES;
}

/*************************************************************************************************************************/
/***************************************************** Poubelle **********************************************************/
/*************************************************************************************************************************/

// size_t FastXFile::getSeqNb() { //Segmentation fault...
//   size_t nbSeq = 0; // Le nombre de seq dans le fichier... Ce que l'on cherche à obtenir
//   char buffer[256];
//   ifs.read(buffer, 256);
//   size_t nb = ifs.gcount();
//   cerr << "nb = " << nb << endl; //DEBUG
//   size_t i = 0;
//   //cerr << "caractère actuel : '" << buffer[i] << "' caractère à venir : '" << buffer[i+1] << "'." << endl; //DEBUG
//   char actual_char = buffer[i];
//   char previous_char = '\n';
//   if ((actual_char == '>' || actual_char == ';' || actual_char == '@') && previous_char == '\n') { //Si le premier caractère du fichier correspond a un début de séquence.
//     ++nbSeq; //On commence a compter les séquences.
//     ++i;
//     cerr << "nbSeq = " << nbSeq << "." << endl;
//     cerr << "i = " << i << "." << endl;
//     cerr << "caractère actuel : '" << actual_char << "' caractère précédent : '" << previous_char << "'." << endl;
//   }
//   while (i < nb) {
//     while(actual_char != '>' || actual_char != ';' || actual_char != '@') {
//       ++i;
//       previous_char = actual_char;
//       actual_char = buffer[i];
//       cerr << "nbSeq = " << nbSeq << "." << endl;
//       cerr << "i = " << i << "." << endl;
//       cerr << "caractère actuel : '" << actual_char << "' caractère précédent : '" << previous_char << "'." << endl;
//     }
//     if ((actual_char == '>' || actual_char == ';' || actual_char == '@') && (previous_char == '\n')) {
//       ++nbSeq;
//       cerr << "nbSeq = " << nbSeq << "." << endl;
//       cerr << "i = " << i << "." << endl;
//       cerr << "caractère actuel : '" << actual_char << "' caractère précédent : '" << previous_char << "'." << endl;
//     }
//     ++i;
//     previous_char = actual_char;
//     actual_char = buffer[i];
//   }
//   return nbSeq;
// }

// FastXSeq FastXFile::parse() {
//   size_t* entetePos(0); 
//   size_t* seqPos(0);
//   size_t* qualityPos(0);
//   char buffer[256];
//   char previous_char = '\n';
//   bool new_sequence = true;
//   bool in_header = true;
//   bool ok = true;
//   size_t nb_nucl = 0;
//   while (ifs && ok) {
//     ifs.read(buffer, 256);
//     size_t nb = ifs.gcount();
//     size_t i = 0;
//     if (nb) {
//       while(buffer[i] < ' ') { //virer les caractères vide en début de fichier.
//         ++i;
//       }
//       //      clog << "Traitement du buffer = '" << buffer << "'" << endl;
//       if (format == FASTA){ //Traitement des fichiers FASTA
//         if (new_sequence) {
//           assert(buffer[i] == '>' || buffer[i] == ';'); //(buffer[0] == '>' || buffer[0] == ';' || buffer[0] == '@') && previous_char == '\n'); 
//           in_header = true;
//           entetePos = ifs.tellg(); // Enregistrer l'adresse dans un pointeur pour entete. peut être même vers le caractère suivant pour éliminer le '>'
//         }
//         while (in_header && (i < nb)) {
//           previous_char = buffer[i];
//           in_header = (buffer[i++] != '\n');
//         }
//         if (!in_header) {
//           new_sequence = false;
//           while (ok && (i < nb)) {
//             if (buffer[i] < ' ') {
//               // Je ne fais rien, il s'agit d'un caractère vide (retour charriot, tab, etc.)
//             } 
//             else {
//               if ((buffer[i] == '>') || (buffer[i] == ';')) { //Ce sont les caractères qui indiquent le début d'une nouvelle séquence
//                 if (previous_char == '\n') {
//                   // il faut que je rembobine ifs de (256 - i) positions.
//                   if (ifs.good()) {
//                     ifs.seekg(-(nb - i), ios_base::cur);
//                   } else {
//                     ifs.clear();
//                     ifs.seekg(-(nb - i), ios_base::end);
//                   }
//                   ok = false;
//                 } else {
//                   // erreur
//                   cerr << "Warning: le caractère '"
//                       << buffer[i]
//                       << "' ne devrait pas apparaître ici"
//                       << endl;
//                   throw "Caractere invalide";
//                 }
//               } 
//               else {
//                 if (Common::isNucl(buffer[i]) || Common::isAA(buffer[i])) { //Deux fonctions qui verifie que les caractères contenus dans les séquences sont bien soit des Acides Aminés, soit des Nucléotides
//                   //Il faut que je stocke seqPos, adresse du premier nucleotide de la sequence... Peut-être il serait plus malin de stocker la dernière position puisqu'on retrouve la première avec le nombre de nucl : nb_nucl.
//                   ++nb_nucl;
//                   // un traitement quelconque, encodage de la seq de ref ?
//                 } else {
//                   // erreur
//                   cerr << "Warning: le caractère '"
//                       << buffer[i]
//                       << "' ne devrait pas apparaître ici"
//                       << endl;
//                   throw "Caractere invalide";
//                 }
//               }
//             }
//             previous_char = buffer[i++];
//           }
//         }
//       }
//       else if (format == FASTQ){ //Traitement des fichiers FASTQ
//         if (new_sequence) {
//           assert(buffer[0] == '@' && previous_char == '\n'); //début de l'entete si @ en début de ligne (problème avec @ c'est si une quality débute par @)
//           in_header = true;
//           //cerr << "entete FASTQ..." << endl; //DEBUG
//         }
//         while (in_header && (i < nb)) {
//           previous_char = buffer[i];
//           in_header = (buffer[i++] != '\n'); //On teste toute la ligne de header. On fait l'assomption que le header n'est que sur une ligne.
//         }
//         if (!in_header) {
//           // cerr << "Debut de seq FASTQ..." << endl;
//           // cerr << "Maintenant, i est égal à : " << i << "." << endl;
//           // cerr << "actual char est : " << buffer[i] << endl;
//           // cerr << "Boolean : in_header : '" << in_header << "', in_quality : '" << in_quality << "', new sequence : '" << new_sequence << "'." << endl;
//           new_sequence = false;
//           while (ok && (i < nb)) {
//             if (buffer[i] < ' ') { 
//               // Caractères vides... On fait rien.
//             } 
//             else {
//               if (buffer[i] == '@') {
//                 if (previous_char == '\n') {
//                   cerr << "Caractère actuel : " << ifs.tellg() << endl;
//                   // in_quality = false;
//                   //cerr << "un @ de début d'entete a été trouvé..." << endl;
//                   // il faut que je rembobine ifs de (256 - i) positions.
//                   if (ifs.good()) {
//                     ifs.seekg(-(nb - i), ios_base::cur);
//                   } else {
//                     ifs.clear();
//                     ifs.seekg(-(nb - i), ios_base::end);
//                   }
//                   ok = false;
//                 } 
//               } 
//               if ((Common::isNucl(buffer[i]) || Common::isAA(buffer[i])) ){// && !in_quality) {
//                 ++nb_nucl;
//                 // un traitement quelconque
//               }
//               if (buffer[i] == '+' && previous_char == '\n'){ //pb si la qualité commence par un +, on peut éventuellement regarder si buffer[++i] est un \n mais pb si qualité fini par un +
//                 // il faut que j'avance de nb_nucl positions.
//                 i += nb_nucl;
//               }
//               // else {
//               //     // erreur
//               //     cerr << "Warning: le caractère '"
//               //         << buffer[i]
//               //         << "' ne devrait pas apparaître ici"
//               //         << endl;
//               //     throw "Caractere invalide";
//             }
//             previous_char = buffer[i++];
//           }
//         }
//       }  
//       else { cerr << "Unsupported format." << endl; // Peut-être facultatif si contrôle avant dans le main.
//       }
//     }
//   }
//   FastXSeq sequenceCurrent(entetePos, seqPos, qualityPos, nb_nucl); //qualityPos reste à une adresse 0 dans le cas du fichier Fasta (pas de quality dans les Fasta).
//   return sequenceCurrent;
// }

#ifdef TEST_FASTXFILE_MAIN
int main(int argc, char** argv) {
  if (argc != 2) {
    cerr << "usage: " << argv[0] << " <fichier>" << endl;
    return 1;
  }

  FastXFile fxf(argv[1]);

  cout << "Format "
       << fxf.getFormat()
       << " détecté." << endl;

  if (fxf.getFormat() == FastXFile::FASTA) {
    size_t nb = 0;
    while (fxf.hasSequence()) {
      cerr << "l'entete est '" << fxf.getCurrentSequenceHeader() << "'" << endl;
      cerr << "la séquence courante est de longueur " << fxf.processCurrentSequence() << endl;
      ++nb;
    }
    cerr << "Il y a " << nb << " sequences" << endl;
  }
  
  cerr << "That's All, Folks!!!" << endl;
  return 0;

}
#endif

#include "fastxseq.h"
#include "Common.h"
#include <string>
#include <iostream>
#include <fstream>
#include <ios>
using namespace std;

//constructeurs

FastXSeq::FastXSeq():entete(NULL), sequence(NULL), quality(NULL), nb_nucl(0){} //vide
FastXSeq::FastXSeq(streampos entete, streampos sequence, streampos quality, size_t nb_nucl):entete(entete), sequence(sequence), quality(quality), nb_nucl(nb_nucl){} //FASTQ
FastXSeq::FastXSeq(streampos entete, streampos sequence, size_t nb_nucl):entete(entete), sequence(sequence), quality(NULL), nb_nucl(nb_nucl){} //FASTA

FastXSeq::~FastXSeq(){
    //delete...
} //implementation du destructeur

bool FastXSeq::parseA(ifstream ifs, streampos entete, streampos sequence) { //bool parseQ(std::streampos entete, std::streampos sequence) { //
  
  this -> entete = entete;
  this -> sequence = sequence;

  //Variables nécessaires au parsing :
  size_t nb_nucl = 0;
  char buffer[4096]; //attention si l'entete est trop long, plus de la taille du buffer, ça chie
  char previous_char = '\n';
  bool ok = true;

  //je me met au début de la seq
  ifs.seekg(sequence); 
  
  while (ifs && ok) {
    ifs.read(buffer, 4096);
    size_t nb = ifs.gcount();
    size_t i = 0;
    while (ok && (i < nb)) {
      if (buffer[i] < '!') {
        // Je ne fais rien, il s'agit d'un caractère vide (retour charriot, tab, etc.)
      } 
      else {
        if ((buffer[i] == '>') || (buffer[i] == ';')) {
          if (previous_char == '\n') {
            // il faut que je rembobine ifs de (4096 - i) positions.
            if (ifs.good()) {
              ifs.seekg(-(nb - i), ios_base::cur);
            } else {
              ifs.clear();
              ifs.seekg(-(nb - i), ios_base::end);
            }
              ok = false;
            } else {
              // erreur caractère invalide.
              cerr << "Warning: character '"
                   << buffer[i]
                   << "' is not a IUPAC character."
                   << endl;
                throw "Invalid character.";
          }
        } 
        else {
          if (Common::isNucl(buffer[i]) || Common::isAA(buffer[i])) {
            ++nb_nucl;
            //Si(bool RefSeq)Alors encode dans un objet type vector<uint8_t>
            //un traitement quelconque, encodage de la seq de ref ?
          } else {
            // erreur caractère invalide.
            cerr << "Warning: character '"
                 << buffer[i]
                 << "' is not a IUPAC character."
                 << endl;
              throw "Invalid character.";
          }
        }
      }
      previous_char = buffer[i++];
    }
  }
  this -> nb_nucl = nb_nucl;
  return !ok;
}

bool FastXSeq::parseQ(ifstream ifs, streampos entete, streampos sequence) { // bool FastXSeq::parseQ(streampos entete, streampos sequence){ //
  
  this -> entete = entete;
  this -> sequence = sequence;

  //Déclaration du streampos pour la position de la qualité
  streampos qualityPos(NULL);

  //Variables nécessaires au parsing :
  size_t nb_nucl = 0;
  char buffer[4096]; //attention si l'entete est trop long, plus de la taille du buffer, ça chie
  char previous_char = '\n';
  bool ok = true;
  bool in_quality = false;
  bool new_sequence = true;
  ifs.seekg(sequence); //je me met au début de la seq

  //je me met au début de la seq
  ifs.seekg(sequence); 
  
  while (ifs && ok) {
    ifs.read(buffer, 4096);
    size_t nb = ifs.gcount();
    size_t i = 0;
    while (ok && (i < nb)) {
      if (buffer[i] < '!') { 
      // Caractères vides... On fait rien.
      } 
      else {
        if (buffer[i] == '+' && previous_char == '\n') {
          //+1 ou +2 pour aller directement sur le premier caractère de quality ?
          qualityPos = (streamoff)ifs.tellg() - (nb - i); 
          if (ifs.good()) {
            ifs.seekg(-(nb - i), ios_base::cur);
          } else {
            ifs.clear();
            ifs.seekg(-(nb - i), ios_base::end);
          }
          in_quality = true;
          new_sequence = false;
        } 
        if ((Common::isNucl(buffer[i]) || Common::isAA(buffer[i])) && !in_quality) { 
          ++nb_nucl;
          // un traitement quelconque
        }
        if (buffer[i] == '@' && previous_char == '\n'){ 
          if (ifs.good()) {
            ifs.seekg(-(nb - i), ios_base::cur);
          } else {
            ifs.clear();
            ifs.seekg(-(nb - i), ios_base::end);
          }
            in_quality = false;
            new_sequence = true;
            ok = false;
          }
          else {
            // erreur caractère invalide.
            cerr << "Warning: character '"
                 << buffer[i]
                 << "' is not a IUPAC character."
                 << endl;
            throw "Invalid character.";
          }
        }
        previous_char = buffer[i++];
      }
    }
  this -> nb_nucl = nb_nucl;
  this -> quality = qualityPos;
  return !ok;
}

string FastXSeq::Seq2String(ifstream ifs){
  char stringSeq[nb_nucl];
  ifs.seekg(entete);
  ifs.read(stringSeq, nb_nucl);
  return stringSeq;
}

ostream &operator<<(ostream &os, FastXSeq fastXSeq) {
  return os << fastXSeq.Seq2String(ifs);
}
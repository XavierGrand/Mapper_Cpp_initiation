#ifndef __FASTXSEQ_H__
#define __FASTXSEQ_H__

#include "EncodedSequence.h"
#include <fstream>
#include <iostream>
#include <ios>
#include <string>

class FastXSeq {

/****************** Partie attributs ****************************************/

 private:
  std::streampos entete; //Je veux stocker l'adresse de l'entete pour faire un getline et recupérer l'entete en string.
  std::streampos sequence; //Je veux stocker l'adresse de la séquence pour pouvoir la récupérer avec la longueur (nb_nucl).
  std::streampos quality; //Je veux stocker l'adresse de la quality, si elle existe sinon initialisée à 0, pour pouvoir la récupérer avec la longueur qui est la même que sequence.
  size_t nb_nucl; //longueur de la séquence qui est obtenue avec la fonction processCurrentSequence de la classe FastXFile

/************************************ Partie Méthodes ***********************/

 public:

  FastXSeq(); //constructeur vide
  FastXSeq(std::streampos entete, std::streampos sequence, std::streampos quality, size_t nb_nucl); //constructeur de FASTQ avec les adresses des attributs d'une seq au format FASTQ
  FastXSeq(std::streampos entete, std::streampos sequence, size_t nb_nucl); //constructeur de FASTA avec les adresses des attributs d'une seq au format FASTA
  ~FastXSeq(); //Si on fait des appels dynamique il faut conserver le destructeur

  /* On renvoie le bool Ok = true si le parsage a réussi */
  bool parseA(std::ifstream ifs, std::streampos entete, std::streampos sequence); // bool parseA(std::streampos entete, std::streampos sequence); //
  bool parseQ(std::ifstream ifs, std::streampos entete, std::streampos sequence); // bool parseQ(std::streampos entete, std::streampos sequence); //

  // Méthode qui récupère la séquence en string.
  std::string Seq2String(ifstream ifs);

  // Getters :
  inline size_t getNbNucl() const {
    return nb_nucl;
  }

  inline streampos getSeqPos() const {
    return sequence;
  }

  inline streampos getQualPos() const {
    return quality;
  }

  inline streampos getEntetePos() const {
    return entete;
  }

};

std::ostream &operator<<(ostream &os, FastXSeq fastXSeq);

#endif

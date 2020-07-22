#ifndef __FASTXFILE_H__
#define __FASTXFILE_H__

#include "fastxseq.h"
#include "EncodedSequence.h"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

class FastXFile {
 
/****************** Partie attributs ******************************************/

 public:

  enum BiologicalSequenceFormat {
    UNKNOWN,
    FASTA,
    FASTQ
  };

 private:
  const char *fname;
  BiologicalSequenceFormat format;
  std::ifstream ifs;
  std::vector<FastXSeq> seqList; //(30) ou (100000)
  /* On pourrait initialiser le vecteur avec 30 positions si fasta de seq de ref, et 100.000 positions si c'est un FastQ de reads.
   *Problème : j'ai pas réussi a faire un shrinkToFit, c++11 nécessaire... */

/********************************** Partie Constructeur ***********************/

 public:

  FastXFile(const char *fname = NULL); //constructeur vide

  ~FastXFile(); //destructeur : pour libérer la mémoire dynamique utilisée lors de la création de l'objet

/************************************ Partie Accesseurs ***********************/

 private:

  void detectFormat(); //setter de format

 public:

  static const char* BioSeqFmt2String(BiologicalSequenceFormat fmt); //getter de BioSeqFormat to string
  
  inline BiologicalSequenceFormat getFormat() const { //getter de format
    return format;
  }

  inline const char *getFilename() const { //getter de filename
    return fname;
  }
  
  std::string getCurrentSequenceHeader(); //getter de Header

/****************************** Partie méthodes *******************************/

  void open(const char* fname); //ouvre le fichier

  void processCurrentSequence(); // FastXSeq //fonction qui parse le fichier pour identifier les différentes séquences

  inline bool hasSequence() const { //tant qu'il reste de la seq à parser...
    return ifs;
  }

  //Méthode qui crée un ojet encodedSeq à partir de la seq lue
  vector<EncodedSequence> EncodedSequence() ; 

/****************************** Partie poubelle *******************************/

  // size_t getSeqNb(); //fonction de comptage des seq, non fonctionnel

  //FastXSeq FastXFile::parse();

};

std::ostream &operator<<(std::ostream &os, FastXFile::BiologicalSequenceFormat fmt); //opérateur <<

#endif

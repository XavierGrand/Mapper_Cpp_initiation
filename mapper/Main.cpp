#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "fastxfile.h"
#include "fastxseq.h"
#include "EncodedSequence.h"
#include "Common.h"
#include "SA_ES.h"
using namespace std;

int main(int argc, char **argv) { // * désigne un pointeur, ** on pointe vers un pointeur d'un array, argc = nb d'arguments, argv = pointer vers les noms des arguments

  /*Check the Usage, ref and reads files*/

  if (argc < 3) {
    //Je prends le parti de dire que la séquence de référence est stockée dans un fasta, classiquement le cas, et les reads dans un fastq.
    cerr << "usage: " << argv[0] << " <refSeq.fasta> <reads.fastq>" << endl; 
    return 1;
  }

  /*Loading Sequences*/

  FastXFile refSeq(argv[1]);
  FastXFile reads(argv[2]);

  size_t nbRefSeq;
  size_t nbReads;

  /*Check the file format, FASTA, FASTQ or UNKNOWN*/

  cout << "Reference Sequence Format : "
       << refSeq.getFormat()
       << " detected." << endl;
  cout << "Reads Sequence Format : "
       << reads.getFormat()
       << " detected." << endl;

  if (refSeq.getFormat() == FastXFile::FASTA) {
    size_t nb = 0;
    while (refSeq.hasSequence()) {
      cerr << "Header : '" << refSeq.getCurrentSequenceHeader() << "'" << endl;
      refSeq.processCurrentSequence();
      ++nbRefSeq;
    }
    cerr << "Reference sequences in file : " << nbRefSeq << " sequences" << endl;
  }
  else { cerr << "usage: " << argv[0] << " <refSeq.fasta> <reads.fastq>" << endl;
  }
  
  if (reads.getFormat() == FastXFile::FASTQ) {
    size_t nb = 0;
    while (reads.hasSequence()) {
      cerr << "Header : '" << reads.getCurrentSequenceHeader() << "'" << endl;
      reads.processCurrentSequence();
      ++nbReads;
    }
    cerr << "Reads in file : " << nbReads << " sequences" << endl;
  }
  else { cerr << "usage: " << argv[0] << " " << argv[1] << " <reads.fastq>" << endl;
  }

  // Le parsing du fichier FASTA, sequence de ref, à produit un vecteur de FastXSeq
  // Passer tous les éléments un par un et les encoder.

  

  cerr << "Mais le code c'est pas 'Le Code' ?" << endl;
  return 0;

} // End of main


  /**************************************************************************************************************************/
  /**************************************************************************************************************************/
  /**************************************************************************************************************************/

  /* Essaie de faire marcher la classe FastXSeq*/

  // string header("SequenceID:1");
  // string seq("ATCG");
  // size_t nb_nuc = 100;
  // cout << header << seq << " nb nucl = " << nb_nuc << endl;
  // FastXSeq testFastXSeq(header, seq, nb_nuc);
  // string test = testFastXSeq.getPosSeq();
  // cout << "J'essaye d'afficher la sequence à partir de l'objet FastXSeq : '" << test << "'." << endl;

  /**************************************************************************************************************************/
  /**************************************************************************************************************************/
  /**************************************************************************************************************************/

  // struct RefSeq {
  //   string entete;
  //   streampos addressRefSeq;
  // };

  // struct ReadSeq {
  //   streampos entete;
  //   streampos addressRead;
  //   streampos quality = NULL;
  // };

//   vector<size_t> addressRefSeq(30); //initialisation du vecteur avec 30 cases mémoires (30 chromosomes)
//   vector<size_t> addressEnteteRefSeq(30); //initialisation du vecteur avec 30 cases mémoires (30 chromosomes)

//   if (refSeq.getFormat() == FastXFile::FASTA) {
//     size_t nb = 0;
//     while (refSeq.hasSequence()) {
//       cerr << "Header : '" << refSeq.getCurrentSequenceHeader() << "'" << endl;
//       cerr << "Sequence length : " << refSeq.processCurrentSequence() << endl;
//       ++nb;
//     }
//     cerr << "Reference sequences in file : " << nb << " sequences" << endl;
//   }
//   else { cerr << "usage: " << argv[0] << " <refSeq.fasta> <reads.fastq>" << endl;
//   }

//   //addressSeq.shrink_to_fit(); //error: ‘class std::vector<long unsigned int>’ has no member named ‘shrink_to_fit’ c++11 needed
    
//   for (size_t i = 0; i < addressRefSeq.size(); i++) {
//     cout << "L'adresse de la séquence " << i + 1 << " est : '" << addressRefSeq[i] << "'." << endl;
//   }

// //   //Encodage

// //   //Suffixes table creation


  
// //   /*Reads Sequences loading*/

//   FastXFile reads(argv[2]);

//   vector<size_t> addressReads; //(100000) reads / Fastq... Dans un premier temps.
//   vector<size_t> addressEnteteReads; //(100000) reads / Fastq... Dans un premier temps.

//   cout << "Reads file Format : "
//        << reads.getFormat()
//        << " detected." << endl;

//   if (reads.getFormat() == FastXFile::FASTQ) {
//     size_t nb = 0;
//     while (reads.hasSequence()) {
//       cerr << "Header : '" << reads.getCurrentSequenceHeader() << "'" << endl;
//       cerr << "Sequence length : " << reads.processCurrentSequence() << endl;
//       ++nb;
//     }
//     cerr << "Reads in file : " << nb << " sequences" << endl;
//   }
//   else { cerr << "usage: " << argv[0] << " " << argv[1] << " <reads.fastq>" << endl;
//   }

//   //addressReads.shrink_to_fit(); //error: ‘class std::vector<long unsigned int>’ has no member named ‘shrink_to_fit’ c++11 needed

//   for (size_t i = 0; i < addressReads.size(); i++) {
//     cout << "L'adresse de la séquence " << i + 1 << " est : '" << addressReads[i] << "'." << endl;
//   }

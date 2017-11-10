#include <zlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

extern "C" {
#include "kseq.h"
}

// convert A to G for hyperediting detection 
// use kseq.h header from heng li 
// http://attractivechaos.github.io/klib/#Kseq%3A%20stream%20buffer%20and%20FASTA%2FQ%20parser

// STEP 1: declare the type of file handler and the read() function
KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[])
{
  gzFile fp;
  kseq_t *seq;

  if (argc == 1) {
    fprintf(stderr, "Usage: %s <original_nt> <replacement_nt> <in.fastq/a.gz> \n", argv[0]);
    return 1;
  }

  fp = gzopen(argv[3], "r"); 
  seq = kseq_init(fp); 
  char ref = *argv[1]; 
  char alt = *argv[2];
  
  char lower_ref = tolower(ref) ; 
  char lower_alt = tolower(alt) ; 
  
  int l ; // error codes are < 0  
  while ((l = kseq_read(seq)) >= 0) { 

    std::string nts = seq->seq.s ; // convert to cpp string
    std::replace(nts.begin(), nts.end(), ref, alt) ;
    std::replace(nts.begin(), nts.end(), lower_ref, lower_alt) ;
    
    if (seq->qual.l == 0){  // FASTA
        std::cout << ">" << seq->name.s << std::endl ;
        std::cout << nts << std::endl ;
    } else if (seq->comment.l > 0) { //FASTQ no index seq
        std::cout << "@" << seq->name.s << " " << seq->comment.s << std::endl ;
        std::cout << nts << std::endl ;
        std::cout << "+" << std::endl ;
        std::cout << seq->qual.s << std::endl ;
    } else { // normal fastq
        std::cout << "@" << seq->name.s << std::endl ;
        std::cout << nts << std::endl ;
        std::cout << "+" << std::endl ;
        std::cout << seq->qual.s << std::endl ;
    }
  }

  kseq_destroy(seq); // STEP 5: destroy seq
  gzclose(fp); // STEP 6: close the file handle

  return 0;
}


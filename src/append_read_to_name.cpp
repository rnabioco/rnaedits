#include <zlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>

extern "C" {
#include "kseq.h"
}

// append read seq to fastq header, so that it can be extracted from bam
// as needed for hyperediting pipeline
//
// use kseq.h header from heng li 
// http://attractivechaos.github.io/klib/#Kseq%3A%20stream%20buffer%20and%20FASTA%2FQ%20parser

// STEP 1: declare the type of file handler and the read() function
KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[])
{
  gzFile fp;
  kseq_t *seq;

  if (argc == 1) {
    fprintf(stderr, "Usage: %s <in.fastq/a.gz> \n", argv[0]);
    return 1;
  }

  fp = gzopen(argv[1], "r"); 
  seq = kseq_init(fp); 

  int l ; // error codes are < 0  
  while ((l = kseq_read(seq)) >= 0) { 

    std::cout << "@" << seq->name.s << ":" << seq->seq.s << " " << seq->comment.s << std::endl ;
    std::cout << seq->seq.s << std::endl ;
    std::cout << "+" << std::endl ;
    std::cout << seq->qual.s << std::endl ;
    
  }

  kseq_destroy(seq); // STEP 5: destroy seq
  gzclose(fp); // STEP 6: close the file handle

  return 0;
}


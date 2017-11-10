#!/usr/bin/env python3
import argparse
import os, sys
import pysam
import gzip
from collections import Counter
from multiprocessing import Pool
from functools import partial

"""
For the supplied intervals, compute the position in 
each read that overlaps the interval. Used to check if
identified editing sites are enriched in a particular 
position across the reads. If enriched at the 5' end
then could be called an edit due to random hexamer mispriming.

See:
Comment on Widespread RNA and DNA Sequence Differences in the Human Transcriptome
Wei Lin1,*, Robert Piskol2,*, Meng How Tan2, Jin Billy Li2, Science
2012.    
http://science.sciencemag.org/content/335/6074/1302.5.full    
"""

def get_edit_positions(bed, bam, qual_threshold, report_mismatch):
    
    bamfile_obj = pysam.AlignmentFile(bam, "rb" )

    pos_dict = Counter() 
    
    for rec in bed:
        contig = rec.split("\t")[0]
        if contig.startswith("#"):
            continue
        start = int(rec.split("\t")[1])
        end = int(rec.split("\t")[2])
        ref = rec.split("\t")[3]
        
        pos_dict = count_reads(bamfile_obj, 
                contig, 
                start,
                end,
                pos_dict,
                qual_threshold,
                report_mismatch)
       

    for pos in sorted(pos_dict):
        print(pos, pos_dict[pos], sep = "\t")

    bamfile_obj.close()

def count_reads(bam, contig, start, end, base_dict, qual_threshold,
        report_mismatch = True):
    
    for pileup in bam.pileup(contig, start, end, max_depth = 100000,
            stepper = "all", truncate = True):

      for read in pileup.pileups:
          # if read has a deletion or is part of N in CIGAR, none is
          # reported for query position
          # ignore multimappers
          NH = read.alignment.get_tag("NH")
          if NH > 1:
              continue
          
          if not read.is_del and not read.is_refskip:
              
              qualities = read.alignment.query_qualities
              
              if qualities[read.query_position] < qual_threshold:
                  continue
              
              if report_mismatch:
                
                # returns list of query position, genomic position, and ref
                # nucleotide
                matches = read.alignment.get_aligned_pairs(matches_only = True, with_seq = True)
                
                match_info = [x for x in matches if x[0] == read.query_position]
                
                # if mismatched then the nucleotide is lowercase
                if match_info[0][2].islower():
                    base_dict[read.query_position] += 1 
                else:    
                    continue

              else:
                base_dict[read.query_position] += 1 

    return base_dict

def main():
    
    parser = argparse.ArgumentParser(description="""This scripts
    For the supplied intervals, compute the position in 
    each read that overlaps the interval. Used to check if
    identified editing sites are enriched in a particular 
    position across the reads. If enriched at the 5' end
    then could be called an edit due to random hexamer mispriming.

    See:
    Comment on Widespread RNA and DNA Sequence Differences in the Human Transcriptome
    Wei Lin1,*, Robert Piskol2,*, Meng How Tan2, Jin Billy Li2, Science
    2012.    
    http://science.sciencemag.org/content/335/6074/1302.5.full""") 
                                      
    parser.add_argument('-b',
                          '--bam',
                          help ='indexed bam file input',
                       required = True)
    parser.add_argument('-v',
                          '--input_bed',
                          help ='input bed file)',
                       required = True)

    parser.add_argument('-q',
                        '--q_threshold',
                        help = """minimum base quality score required for
                        counting base in overlapping read""", 
                        default = 20, type = int)

    parser.add_argument('-m',
                        '--report_only_mismatches',
                        help = """report positions of mismatches only,
                        rather than just position of site within read""", 
                        action = 'store_true')

    args=parser.parse_args()
    
    gzopen = lambda f: gzip.open(f, 'rt') if f.endswith(".gz") else open(f)
    
    bam_name = args.bam
    bed = gzopen(args.input_bed)
    
    print("position\ttotal_reads")
    
    get_edit_positions(bed, 
            bam_name, 
            qual_threshold = args.q_threshold,
            report_mismatch = args.report_only_mismatches)

if __name__ == '__main__': main()



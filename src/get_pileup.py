#!/usr/bin/env python3
import argparse
import os, sys
import pysam
import gzip
from collections import Counter
from multiprocessing import Pool
from functools import partial

"""
Parse vcf and count number of reference and alternative alleles for each
variant. Requires indexed bam file. The alternative alleles reported will
only be SNPs. 
"""

def base_counter():
    """create counter object populated with canonical bases"""
    
    bases_to_count = ["A","T","C","G"]
    base_dict = Counter()
    
    for base in bases_to_count:
        base_dict[base] = 0  
    
    return base_dict

def get_allele_counts(vcf, bam, qual_threshold, 
        trimmed_position = 0, multithread = False):
     
    # build list of results if multithreaded, otherwise just print
    if multithread:
        out_results = []

    bamfile_obj = pysam.AlignmentFile(bam, "rb" )
    for rec in vcf:
        contig = rec.split("\t")[0]
        if contig.startswith("#"):
            continue
        pos = int(rec.split("\t")[1])
        ref = rec.split("\t")[3]
    
        # only keep snps
        if len(ref) > 1:
            continue

        base_dict = base_counter() 
        base_dict = count_reads(bamfile_obj, 
                contig, 
                pos, 
                base_dict,
                qual_threshold,
                trimmed_position)
        
        if multithread:
            result_builder(contig, pos, ref, base_dict, out_results)
        else:
            count_writer(contig, pos, ref, base_dict)

        base_dict.clear()
    
    if multithread: 
        return out_results

    bamfile_obj.close()

def count_reads(bam, contig, pos, base_dict, qual_threshold,
        trimmed_position = 0):
    """count the number of reads a single base. 

    Reads will not be counted if they are not mapped, secondary
    alignments, failed QC, marked duplicate, set as multiply mapped (NH >
    1), or if the overlapping base quality does not pass a
    quality_threshold. 
    
    The coverage is reported per base for A, T, C, and G. N in the read
    are ignored.

    Parameters
    ----------

    bam : pysam.AlignmentFile object
        Open file object with matching bam index file

    contig: string
        chromosome/contig name 

    pos: int
        one-based position of the base to be counted

    base_dict: dictionary
        dictionary of bases desired to be counted, initialized with keys
        of bases to be counted and zero counts

    qual_threshold: int
        minimum base quality required for counting base
    
    trimmed_position: int
        exclude 5' nucleotides up to this position, default = 0 (no
        trimming)

    Returns
    -------

    a dictionary with keys A,T,C,G with values of counts

    """

    for pileup in bam.pileup(contig, pos - 1, pos, max_depth = 100000, 
            stepper = "all", truncate = True):

      for read in pileup.pileups:
          # if read has a deletion or is part of N in CIGAR, none is
          # reported for query position
          # ignore multimappers
          NH = read.alignment.get_tag("NH")
          if NH > 1:
              continue
          
          if not read.is_del and not read.is_refskip:
              
              # do not count alleles at 5' end of read 
              # may be due to random hexamer mispriming

              if read.query_position < trimmed_position:
                  continue

              qualities = read.alignment.query_qualities
              
              if qualities[read.query_position] < qual_threshold:
                  continue
              
              found_base = read.alignment.query_sequence[read.query_position]
              
              if found_base in base_dict:
                  base_dict[found_base] += 1
              else:
                  continue

    return base_dict

def count_writer(contig, pos, ref, base_dict):
    for bases, counts in base_dict.items():
        print(contig, pos, ref, bases, counts, sep = "\t")

def result_builder(contig, pos, ref, base_dict, list_obj):
    """ format results into list object for multithreaded"""
    for bases, counts in base_dict.items():
        res = "{}\t{}\t{}\t{}\t{}".format(contig, pos, ref, bases, counts)
        list_obj.append(res)

def emit_chunks(chunk_size, file_path):
    # http://stackoverflow.com/questions/31164731/python-chunking-csv-file-multiproccessing
    lines_count = 0
    gzopen = lambda f: gzip.open(f, 'rt') if f.endswith(".gz") else open(f)
    
    with gzopen(file_path) as f:
        chunk = []
        for line in f:
            lines_count += 1
            chunk.append(line)
            if lines_count == chunk_size:
                lines_count = 0
                yield chunk
                chunk = []
            else:
                continue
        if chunk : yield chunk


def main():
    
    parser = argparse.ArgumentParser(description="""This scripts
    calculates the number of alternative bases observed for each supplied
    region from a vcf file. Requires an indexed bam. Only SNPs are
    reported, and by default reads are not counted if they are secondary,
    duplicates, unmapped, failedQC, or below a base quality threshold""")
                                      
    parser.add_argument('-b',
                          '--bam',
                          help ='indexed bam file input',
                       required = True)
    parser.add_argument('-v',
                          '--input_vcf',
                          help ='input vcf file)',
                       required = True)
    parser.add_argument('-t',
                        '--threads',
                        help = 'number of threads to use, defaults to 1', 
                        default = 1, type = int)

    parser.add_argument('-q',
                        '--q_threshold',
                        help = """minimum base quality score required for
                        counting base in overlapping read""", 
                        default = 20, type = int)
    
    parser.add_argument('-p',
                        '--trim_5p',
                        help = """trim 5p nucleotides up to this point, 
                        default is 0 which is no trimming""", 
                        default = 0, type = int)

    args=parser.parse_args()
    
    gzopen = lambda f: gzip.open(f, 'rt') if f.endswith(".gz") else open(f)
    
    bam_name = args.bam
    
    #print header
    print("#chrom", "pos", "ref", "base", "counts", sep = "\t") 

    if args.threads > 1:
    
        # chunk up vcf file into 10,000 records
        chunk_size = 10000 
        gen = emit_chunks(chunk_size, args.input_vcf)

        # set up processes
        pool = Pool(args.threads)

        # build function with multiple args
        func = partial(get_allele_counts, 
                bam = bam_name, 
                qual_threshold = args.q_threshold,
                trimmed_position = args.trim_5p,
                multithread = True)
        
        # count alleles in parallel, printing result in a threadsafe
        # manner
        for res in pool.imap(func, gen):
            for allele in res:
                print(allele)
    else:
        vcf = gzopen(args.input_vcf)
        get_allele_counts(vcf, bam_name, 
                qual_threshold = args.q_threshold,
                trimmed_position = args.trim_5p)

if __name__ == '__main__': main()



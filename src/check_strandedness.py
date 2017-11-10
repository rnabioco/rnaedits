#!/usr/bin/env python3
import argparse
import os, sys
import pysam
import gzip
from collections import Counter

"""
Parse vcf and extract reads overlapping each variant. Enumerate strandedness 
of overlapping reads. Requires indexed bam file. 
"""

def read_quality_check(read):
    """take pysam alignment object and return bool true if passes
    filters for paired end RF orientation data
    """
    
    if read.is_duplicate or read.is_qcfail or read.is_secondary or read.is_supplementary or read.mate_is_unmapped:
        return False
    
    elif read.is_paired and read.is_proper_pair:
        return True

    else:
        sys.exit("error flag parsing incorrect for read {}".format(read.qname))

def strandedness_check(read):
    """take pysam alignment object and return 0 for positive and 1 for
    negative, assumes RF orientation (Read2 on forward is forward)
    """

    strand = -1

    if read.is_read1:
        if read.is_reverse and not read.mate_is_reverse:
            strand = 0 
        elif read.mate_is_reverse and not read.is_reverse:
            strand = 1
        else:
            sys.exit("error strandparsing incorrect for read {}".format(read.qname))

    else:
        if read.is_reverse and not read.mate_is_reverse:
            strand = 1 
        elif read.mate_is_reverse and not read.is_reverse:
            strand = 0
        else:
            sys.exit("error strandparsing incorrect for read {}".format(read.qname))
       
    return strand 

def get_allele_strandedness(bam, vcf):
    """iterate through bam and check for read quality.
    Enumerate positive and negative reads overlapping each variant"""

    samfile = pysam.AlignmentFile(bam, "rb" )

    for rec in vcf:
        contig = rec.split("\t")[0]
        
        if contig.startswith("#"):
            continue
        
        pos = int(rec.split("\t")[1])
        ref = rec.split("\t")[3]
        
        strand_count = {}
        
        for strand in ["positive", "negative"]:
            strand_count[strand] = 0
        
        for read in samfile.fetch(contig, pos - 1, pos):
            
            # ignore multimappers
            NH = read.get_tag("NH")
            if NH > 1:
                continue
            
            if read_quality_check(read) is False :
                continue
            else:
                strand = strandedness_check(read)
                if strand is 0:
                    strand_count["positive"] += 1
                elif strand is 1: 
                    strand_count["negative"] += 1
                elif strand is -1:
                    sys.exit("strand parsing failure at read {}".format(read.qname))
        
        if strand_count["positive"] > 0 or strand_count["negative"] > 0:
            pos_counts = strand_count["positive"]
            neg_counts = strand_count["negative"]

            prop_pos = float(pos_counts) / (pos_counts + neg_counts) 
            prop_neg = float(neg_counts) / (pos_counts + neg_counts)
        
        else:
            #report zeros for sites with no reads overlapping
            pos_counts, neg_counts, prop_pos, prop_neg = [0, 0, 0, 0]
            
        print(contig, pos, ref, pos_counts, neg_counts, prop_pos,
                    prop_neg,  sep = "\t")
    samfile.close()

def main():
    
    parser = argparse.ArgumentParser(description="""This scripts 
    calculates the strandedness of each proper paired-end alignment over
    a set of sites defined by a vcf. Requires an indexed bam. Only SNPs are
    reported, and by default reads are not counted if they are secondary,
    duplicates, unmapped, failedQc, or non-unique based on NH tag""")
                                      
    parser.add_argument('-b',
                          '--bam',
                          help ='indexed bam file input',
                       required = True)
    parser.add_argument('-v',
                          '--input_vcf',
                          help ='input vcf file)',
                       required = True)
    args=parser.parse_args()
    
    gzopen = lambda f: gzip.open(f, 'rt') if f.endswith(".gz") else open(f)
    
    bam = args.bam
    vcf = gzopen(args.input_vcf)
     
    get_allele_strandedness(bam, vcf)
    vcf.close()

if __name__ == '__main__': main()




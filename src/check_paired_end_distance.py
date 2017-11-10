#!/usr/bin/env python3
import argparse
import os, sys
import pysam
import gzip
from collections import Counter

def process_bam(bam, bedfile, editfile, bed_out, edit_out):
    
    bamfile_obj = pysam.AlignmentFile(bam, "rb" )
    bedfile_obj = open(bedfile, 'r')
    editfile_obj = open(editfile, 'r')

    bed = open(bed_out, 'w')
    edits = open(edit_out, 'w')

    chrom_names = bamfile_obj.references
    chrom_lens = bamfile_obj.lengths
    
    good_reads = []
    for line in bedfile_obj:
        
        interval = line.split("\t")

        chrom, start, end, read = interval[:4]
        start = int(start)
        end = int(end)

        # determine max chrom length
        max_len = chrom_lens[chrom_names.index(chrom)]
        
        # +/- 500kbp or ends of chrom
        
        slop_start = max([0, start - 500000])
        slop_end = min(end + 500000, max_len)
                
        # query bam for interval +/- 500kbp
        # check for read name in this region (probably faster than
        # iterating over entire bam for each record)
        

        for alignment in bamfile_obj.fetch(reference = chrom, 
                                           start = slop_start,
                                           end = slop_end):
            bamread = alignment.query_name
            # strip off last field (sequence)
            bamread = bamread[:bamread.rfind(':')]
            
            if bamread == read:
                good_reads.append(read)
                bed.write(line)

    for line in editfile_obj:
        edit_read = line.split("\t")[3]
        if edit_read in good_reads:
            edits.write(line)


    bamfile_obj.close()
    bedfile_obj.close()
    editfile_obj.close()
    bed.close()
    edits.close()


def main():
    
    parser = argparse.ArgumentParser(description="""
    this script takes a bed file with read position and name
    and searches for the mate in the bam. If mate is aligned, 
    and within 500kbp then the read is returned. 
    """)

    parser.add_argument('-b',
                          '--bam',
                          help ='indexed bam file input',
                       required = True)
    parser.add_argument('-r',
                          '--reads',
                          help ="""bedfile with reads, readname needs to be
                          in fourth column""",
                       required = True)
    parser.add_argument('-e',
                        '--edits',
                          help ="""bedfile with edits, readname needs to be
                          in fourth column""",
                       required = True)
    parser.add_argument('-R',
                          '--reads_out',
                          help ="""output file for reads""",
                       required = True)
    parser.add_argument('-E',
                          '--edits_out',
                          help ="""output file for edits""",
                       required = True)
                        
    args=parser.parse_args()
    
    bam_name = args.bam
    reads_name = args.reads
    edits_name = args.edits
    process_bam(bam_name, reads_name, edits_name, args.reads_out,
            args.edits_out)

if __name__ == '__main__': main()



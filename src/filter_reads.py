#!/usr/bin/env python3
import argparse
import os, sys
import gzip
from collections import Counter
from itertools import islice, groupby
from scipy.stats.mstats import mquantiles

class FastqRecord:
    '''storage class for fastq data see https://github.com/brwnj/umitools/blob/master/umitools/umitools.py'''

    def __init__(self, args):
        self.header = args[0][1:]
        self.seq = args[1]
        self.qual = args[3]
        assert len(self.seq) == len(self.qual), "Seq and Qual vary in length"

    def __repr__(self):
        return "Fastq({})".format(self.header)

    def __str__(self):
        return "@{}\n{}\n+\n{}".format(self.header, self.seq, self.qual)


def parse_fq(file_handle):
    ''' read filehandle and yield fastq record '''
    
    fqclean = (x.strip("\r\n") for x in file_handle if x.strip())
    while True:
        rd = [x for x in islice(fqclean, 4)]
        if not rd:
            raise StopIteration
        assert all(rd) and len(rd) == 4
        yield FastqRecord(rd)

def base_counter():
    """create counter object populated with canonical bases"""
    
    bases_to_count = ["A","T","C","G","N"]
    base_dict = Counter()
    
    for base in bases_to_count:
        base_dict[base] = 0  
    
    return base_dict

def check_nucleotide_content(seq):
    
    """ check sequence for good nucleotide content 
    no single nucleotide > 60% or < 10%
    no N's > 10% of read
    no strech of a single nt longer than 20
    return False if any conditions are not satisfied
    """
    seq_length = len(seq)
    nt_percent = dict() 
    nt_count = base_counter()
    
    for i in seq:
        nt_count[i] += 1
     
    nt_percent = dict() 
    
    for nt, count in nt_count.items():
        nt_percent[nt] = 100 * (count / seq_length)

    for nt, percent in nt_percent.items():
        if nt == "N" and percent > 10:
            return False
        if nt == "N":
            continue
        if percent > 60 or percent < 10:
            return False

    #compute run length encoding for seq
    seq_rle = [(len(list(count)), nt) for nt, count in groupby(seq)]
    
    # if single nucleotide is repeated for > 20 positions drop read
    for unique_seq in seq_rle:
        if unique_seq[0] > 20:
            return False
    return True

def check_qual(qual):
    """check seq for mean quality < 25 
       drop qual scores in lowest decile
    """
    
    quals = []
    for i in qual:
        quals.append(ord(i) - 33)

    # drop lowest decile of qual scores
    decile = float(mquantiles(quals, prob = [0.1]))
    
    quals = [x for x in quals if x > decile]

    mean_qual = float(sum(quals)) / max(len(quals), 1)
    
    if mean_qual < 25:
        return False
    else:
        return True

def check_repeats(seq, repeats, num_repeats = 10):
    
    """ check for repeat but exclude single nucleotide repeats """

    for test_pattern in repeats:     
        if len(test_pattern) is 1:
            continue

        test_pattern_rpt = test_pattern * num_repeats

        if test_pattern_rpt in seq:
            return False

    return True 

def write_fq(read, outfile_handler):
    """ write fastq to disk gzipped """
        
    outfile_handler.write(str(read) + "\n")

def process_fastq(fastq_file, repeats_file, gzip_fh_out):
    
    test_patterns = [x.strip('\r\n') for x in repeats_file]

    for read in parse_fq(fastq_file):
        
        # if single nucleotide is > 60% or < 10% drop read
        # or if N content is > 10%
        # or single nt repeated > 20 in a row
        if not check_nucleotide_content(read.seq):
            continue
        
        # drop low quality reads
        if not check_qual(read.qual):
            continue
        
        # drop reads with high levels of repeats
        if not check_repeats(read.seq, test_patterns):
            continue

        else:
            write_fq(read, gzip_fh_out)

def main():
    
    parser = argparse.ArgumentParser(description="""
    Filter fastq (optionally gzipped) to remove possible
    sequencing artifacts. Used for prefiltering reads
    for detecting hyperedited RNA-Seq reads.
    See: 
    Porath et al. A genome-wide map of hyper-edited RNA reveals numerous
    new sites. Nature Communications. 2014. doi:10.1038/ncomms5726
    """)

    parser.add_argument('-i',
                          '--fastq',
                          help ='fastq to filter',
                       required = True)
    parser.add_argument('-r',
                          '--repeats',
                          help ='repeats_file',
                       required = True)
    parser.add_argument('-o',
                          '--output',
                          help ='output filename for gzipped fastq',
                       required = True)

    args=parser.parse_args()
    
    gzopen = lambda f: gzip.open(f, 'rt') if f.endswith(".gz") else open(f)
    gzwrite = lambda f: gzip.open(f, 'wt') if f.endswith(".gz") else open(f, 'w')    
    
    fq = gzopen(args.fastq)
    rpts = gzopen(args.repeats)
    outf = gzwrite(args.output)
    
    process_fastq(fq, rpts, outf)
    
    outf.close()
    rpts.close()
    fq.close()

if __name__ == '__main__': main()



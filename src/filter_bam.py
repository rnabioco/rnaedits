#!/usr/bin/env python3
import argparse
import os, sys
import pysam
import gzip
from collections import Counter

"""
Parse bams produced to detect hyperedited reads and
filter to enriched for bona fide hyperedited regions
"""


def base_counter():
    """create counter object populated with canonical bases"""
    
    bases_to_count = ["A","T","C","G"]
    base_dict = Counter()
    
    for base in bases_to_count:
        base_dict[base] = 0  
    
    return base_dict

def rev_comp(seq):

    seq = seq.upper()
    
    bps = bp_complement()
    
    comp_seq = ""
    for nt in seq:
        comp_seq += bps[nt]

    return comp_seq[::-1]

def recover_query_seq_and_quals(alignment_obj):
    """ take pysam alignment object and recover original seq
    from query name and return list of phred scores matching seq"""

    qname = alignment_obj.query_name
    qname_elements = qname.split(":")
    qname_seq = qname_elements[-1]
    qquals = alignment_obj.query_qualities
    
    # need to subset to match actual aligned portion
    qstart = alignment_obj.query_alignment_start
    qend = alignment_obj.query_alignment_end
    aligned_seq = qname_seq[qstart:qend]
    aligned_quals = list(qquals[qstart:qend])

    return aligned_seq, aligned_quals

def recover_reference_seq(alignment_obj, fasta_file):
    """ recover original reference sequence from
    fasta file """
    rchrom = alignment_obj.reference_name
    rstart = alignment_obj.reference_start
    rend = alignment_obj.reference_end
    
    r_seq = fasta_file.fetch(reference=rchrom, start= rstart, end= rend)
    return r_seq

def align_seqs(seq1, seq2, map_qual):
    """ simple ungapped alignment between two strings of the same length
    also filter mismatches to only keep high quality mismatches (PHRED
    >=30)
    returns list with query nt and ref nt and position of mismatch
    """

    assert len(seq1) == len(seq2) == len(map_qual)
    
    mismatches = []
    for i in range(len(seq1)):
        if seq1[i] == seq2[i]:
            continue
        else:
            if map_qual[i] >= 30:
                mms = [seq1[i], seq2[i], i]
                mismatches.append(mms)

    return mismatches

def count_edits(alignment_list, read_type, ref, alt, c_ref, c_alt):

    """take output from align_seqs and count a to g
    based on strand return editing count
    also count total mismatches and filter alignment to 
    keep only edits. assumes paired end stranded with RF orientation
    i.e. read 1 is rev. comp to transcribed strand read 2 is sense.
    return both values and filtered alignments"""
    
    edit_count = 0 
    total_count = 0
    
    filtered_alignment = []
    for mm in alignment_list:
        total_count += 1

        q_nt = mm[0]
        r_nt = mm[1]
        if read_type == "R2":
            if q_nt == alt and r_nt == ref:
                edit_count += 1
                filtered_alignment.append(mm)
        else:
            if q_nt == c_alt and r_nt == c_ref:
                edit_count += 1
                filtered_alignment.append(mm)

    return edit_count, total_count, filtered_alignment

def cluster_position(alignment_list):
    """ take alignment list from align_seqs (or count_edits) 
    and find start and end of cluster """
    
    start = alignment_list[0][2]
    end = alignment_list[-1][2]
    assert start < end
    return [start, end]

def filter_cluster_position(start, end, read_length):
    """ return false if start and end within first 20% or
    last 20% of read """
    
    assert start < end

    first_20_pos = round(read_length * 0.20)
    last_20_pos = round(read_length * 0.80)
    
    if end < first_20_pos:
        return False
    elif start > last_20_pos:
        return False
    else:
        return True

def bp_complement():

    nt = {'A':'T',
          'T':'A',
          'G':'C',
          'C':'G',
          'N':'N'}
    return nt

def base_dict():
    
    nt = {'A':0,
          'T':0,
          'G':0,
          'C':0,
          'N':0}

    return nt

def nt_comp(seq):
    
    nt_dict = base_dict()

    for i in seq:
        nt_dict[i] += 1
    
    nt_percent = dict()
    for nt, count in nt_dict.items():
        nt_percent[nt] = 100.0 * (count / len(seq))

    return nt_percent

def best_multimapper(alignment, fasta_obj, read_type, ref, alt, c_ref, c_alt):
    """ iterate through all possible best alignments
    and select alignment with highest % of A to G mismatches
    provided that it is > 10% than next best alignment """
    
    alignment_info = []
    
    # get list of additional alignments
    xa_list = alignment.get_tag("XA")
    xa_list = xa_list.split(";")
    
    # compute percent A to G mismatch over total for
    # reported alignment

    qseq, qquals = recover_query_seq_and_quals(alignment)
    rseq = recover_reference_seq(alignment, fasta_obj)
        
    # get strand of alignment
    if alignment.is_reverse:
        strand = "-"
        rseq = rev_comp(rseq)
    else:
        strand = "+"

    new_align = align_seqs(qseq, rseq, qquals)
    
    edit_count, total_count, filtered_alignment = count_edits(new_align, read_type, ref, alt, c_ref, c_alt)
        
    q_len = alignment.query_alignment_length
    
    percent = edit_count / float(q_len)
    
    alignment_info.append([alignment.reference_name, 
        alignment.reference_start,
        alignment.cigarstring,
        alignment.get_tag("XM"),
        percent])
    
    # iterate through XA defined alignments
    # calculate percent A to G over others

    for other_alignment in xa_list:
        if other_alignment == '': continue 
        chrom, position, cigar, nm = other_alignment.split(",")
        pos = int(position[1:]) #strip off "+"
        strand = position[0] # "+" or "-"

        # parse cigar and get alignment lenght
        # due to bwa mapping settings, this should always be 
        # the lenght of the input read, defined as 125M for ex.
        align_len = int(cigar[:-1])

        r_seq = fasta_obj.fetch(reference=chrom, 
                start = pos, 
                end = pos + align_len)
        if strand == "-": 
            r_seq = rev_comp(r_seq)

        new_align = align_seqs(qseq, rseq, qquals)
        edit_count, total_count, filtered_alignment = count_edits(new_align, read_type, ref, alt, c_ref, c_alt)
        percent = edit_count / float(q_len)
        
        alignment_info.append([chrom, 
          pos,
          cigar,
          nm,
          percent])
    
    # select best alignment if found
    # based on percent a to g
    percent_ag = []
    for align in alignment_info:    
        percent_ag.append(align[4])

    sorted_ag = sorted(percent_ag, reverse = True) 
    
    max_ag = sorted_ag[0]
    second_ag = sorted_ag[1]
    
    if max_ag > (second_ag * 1.10):

        best_alignment_idx = percent_ag.index(max_ag)
        best_alignment = alignment_info[best_alignment_idx]
        alignment.reference_start = best_alignment[1]
        return alignment

def process_bam(bam, fasta, output_edits, output_bed, read_type, edit_type):
    
    bamfile_obj = pysam.AlignmentFile(bam, "rb" )
    fasta_obj = pysam.FastaFile(fasta)

    edits_out = open(output_edits, 'w')
    reads_out = open(output_bed, 'w')

    ref = edit_type[0]
    alt = edit_type[1]
    c_ref = rev_comp(ref)
    c_alt = rev_comp(alt)
    
    for alignment in bamfile_obj.fetch():
         
        # check for uniqueness
        if alignment.has_tag("XA"): 
            best_alignment = best_multimapper(alignment, fasta_obj,
                    read_type, ref, alt, c_ref, c_alt)  
            if best_alignment is not None:
                alignment = best_alignment

        # original query sequence
        qseq, qquals = recover_query_seq_and_quals(alignment)
        
        # reference sequence
        rseq = recover_reference_seq(alignment, fasta_obj)
        
        # get strand of alignment
        if alignment.is_reverse:
            strand = "-"
            rseq = rev_comp(rseq)
        else:
            strand = "+"
        
        # align sequences
        new_align = align_seqs(qseq, rseq, qquals)

        # count editing and total mismatches
        edit_count, total_count, filtered_alignment = count_edits(new_align, read_type, ref, alt, c_ref, c_alt)
        
        q_len = alignment.query_alignment_length
        
        #require at least 5% of sequence is edited
        if (edit_count / float(q_len)) < 0.05:
            continue

        #require 60% of mismatches are a_to_g
        prop_mismatch = (edit_count / float(total_count))
        if prop_mismatch < 0.60:
            continue
        
        #require 80% of mismatches are a_to_g if read length <= than 60
        if q_len <= 60 and prop_mismatch < 0.80:
            continue
        
        cluster_pos = cluster_position(filtered_alignment)
        
        cluster_length = cluster_pos[1] - cluster_pos[0]

        # if cluster is less than 10% of alignment, discard
        if (cluster_length / float(q_len)) < 0.10:
            continue

        # if cluster entirely within first or last 20% of read discard
        good_position = filter_cluster_position(cluster_pos[0], cluster_pos[1], q_len)
        
        if good_position == False:
            continue

        # only keep clusters with < 60% of a single nt
        cluster_seq = qseq[cluster_pos[0]:cluster_pos[1] + 1] # inclusive
    
        cluster_nt_comp = nt_comp(cluster_seq)
         
        for percent in cluster_nt_comp.values():
            if percent > 60:
                continue
        
        # write out bed-like file (zero based half open)
        # one entry per editing site

        chrom = alignment.reference_name
        align_start = alignment.reference_start
        
        read_id = alignment.query_name
        read_id = read_id[:read_id.rfind(":")] #drop seq from end of id
        c_start = align_start + cluster_pos[0]
        c_end = align_start + cluster_pos[1] + 1 #(zero based half open)
        c_pos = "{}:{}".format(c_start, c_end)
       
        # report strand based on Read type (i.e. + strand alignment
        # from R1 is actually derived from RNA on the - strand)
        # R2 reports proper strand, so no modification necessary

        if read_type == "R1":
            if strand == "+":
                strand = "-"
            else:
                strand = "+"

        read_output = [chrom, str(c_start), str(c_end), read_id,
                str(len(filtered_alignment)), strand]

        read_out = "\t".join(read_output)
        reads_out.write("{}\n".format(read_out))

        for mm in filtered_alignment:
            start = align_start + mm[2]
            end = start + 1
            start = str(start)
            end = str(end)
            edit_out = [chrom, start, end, read_id, c_pos, strand]
            edit_out = "\t".join(edit_out)
            edits_out.write("{}\n".format(edit_out))

    bamfile_obj.close()
    edits_out.close()
    reads_out.close()


def main():
    
    parser = argparse.ArgumentParser(description="""
    This script parses a bam file produced for detecting hyperedited reads
    and attempts to remove likely mismapped low quality alignments """ )

    parser.add_argument('-b',
                          '--bam',
                          help ='indexed bam file input',
                       required = True)
    parser.add_argument('-f',
                          '--ref_fasta',
                          help ='reference fasta indexed with samtools faidx',
                       required = True)
    parser.add_argument('-e',
                          '--bededitout',
                          help ='output bed file, one entry per edit',
                       required = True)
    parser.add_argument('-r',
                          '--bedreadout',
                          help ='output bed file, one entry per read',
                       required = True)
    parser.add_argument('-s',
                          '--strandedness',
                          help ="""R1 or R2 from paired end alignment, 
                          determines which mismatch is counted. i.e
                          for R2 in PE RF orientation, an alignment on the
                          reverse strand should contain G's, but R1 should
                          have C's if edited, supply either 'R1' or 'R2' """,
                       required = True)
    parser.add_argument('-a',
                          '--edit_site',
                          help ="""
                          editing type, supply as two letters. i.e.
                          "A" "G" is a to i editing""",
                          nargs="+",
                       required = False, default = ["A", "G"])

    args=parser.parse_args()
    
    bam_name = args.bam
    fa_name = args.ref_fasta
    edits = args.bededitout
    reads = args.bedreadout
    strandedness = args.strandedness
    edit_type = args.edit_site
    process_bam(bam_name, fa_name, edits, reads, strandedness, edit_type)
    
if __name__ == '__main__': main()



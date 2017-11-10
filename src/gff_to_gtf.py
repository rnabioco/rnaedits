#! /usr/bin/env python3

import os
import sys
import argparse
from collections import defaultdict
import tempfile
import re

# from gtf_to_bed.py
from gtf_to_bed import gtf_to_bed

def build_ref_table(gtf, attributes, key_col = 7, region = "transcript"):
    """ return dictionary of transcript ids:list of other metadata """
    
    transcript_map = defaultdict(set)
    
    tmp = tempfile.NamedTemporaryFile(mode = 'w')
    
    # write table to temporary file
    # extract out region field (3) also
    gtf_to_bed(gtf, attributes, [3], tmp)

    # need to rewind file pointer to start
    tmp.seek(0)
    
    nattrs = len(attributes)
    # parse temporary file into transcript map
    for line in open(tmp.name):
        tmp_cols = line.split("\t")
        feature = tmp_cols[-1].rstrip("\n")
        
        # only build map from transcript regions
        if feature != region: 
            continue
        
        #make transcript id key, and build set from attributes
        tid = tmp_cols[key_col - 1]
        metadata_list = [x for x in tmp_cols[7:-1]]
        transcript_map[tid].add(tuple(metadata_list))
    
    tmp.close()

    return transcript_map 

#def reference_lookup(gtf):


def parse_transdecoder_gff(gff, transcript_map):

    for line in gff:
        # ignore header lines, 
        # and blank lines which are present in transdecoder output
        if line.startswith("#") or line.strip() == "":
            continue
        
        fields = line.split("\t")
        # drop the gene regions, use mRNA regions instead
        if fields[2] == "gene":
            continue
        
        # parse out transcript ID 
        attrs = fields[8]
        tid_pattern = re.compile('ID=(cds.Gene|Gene)[.][0-9]+::(TU[0-9]+)::')
        tid = tid_pattern.search(attrs).groups()[1]
        metadata = transcript_map[tid]
        
        if len(metadata) > 1:
            print("more than 1 transcript annotation found for {}".format(tid),
                 file = sys.stderr)

        elif len(metadata) == 0:
            print("""unable to parse out transcript_id and other
                     attributes for line \n{}""".format(line),
                   file = sys.stderr)
        else:
            for attributes in metadata:
                denovo_gene_id = attributes[0]
                ref_gene_id = attributes[1]
        
        outattrs = [denovo_gene_id, tid, ref_gene_id] 
        
        attrs_format = 'gene_id "{d[0]}"; transcript_id "{d[1]}"; gene_name "{d[2]}"; gene_biotype "protein_coding";'
        res = attrs_format.format(d = outattrs)

        # return result, only coding regions are annotated in the
        # transdecoder gff

        print("\t".join(fields[:8]), res, sep = "\t" )

def main():

    parser = argparse.ArgumentParser(description="""
    convert gff3 file produced by Transdecoder into GTF suitable
    for use with SnpEFF. Also add in gene/transcript attributes
    """)
                                      
    parser.add_argument('-i',
                          '--gff3',
                          help ='transdecoder gff3 file',
                       required = True)
    parser.add_argument('-g',
                          '--gtf',
                          help ='original gtf used to produce transdecoder file',
                       required = True)
    parser.add_argument('-a',
            '--attributes', 
            help="""list of attributes in gtf to add to gff. The first
            listed attribute will be used as the lookup key and must match
            the lookup key field taken from the transdecoder gff3 file""",
            required = True, nargs='+')
    args=parser.parse_args()
   
    gff = args.gff3
    gtf = args.gtf
    attributes = args.attributes
    
    # generate lookup table for transcript ids and other metadata from
    # original gtf. Loaded into memory as dictionary
    gtf_obj = open(gtf, 'r')
    gff_obj = open(gff, 'r')

    transcript_map = build_ref_table(gtf_obj, attributes)
    
    parse_transdecoder_gff(gff_obj, transcript_map)

if __name__ == "__main__": main()

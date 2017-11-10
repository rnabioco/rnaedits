#!/usr/bin/env python3
import argparse
import os, sys
import gzip
import re

"""
Parse vcf and extract out data as a BED extended 
"""

def vcf_to_bed(vcf, attrs):
    
    if attrs is not None:
        attrs_len = len(attrs)
    
    for rec in vcf:
        
        rec = rec.split("\t")
        contig = rec[0]
        
        if contig.startswith("#"):
            continue
        
        start = int(rec[1]) - 1
        end = start + 1
        ref = rec[3]
        alt = rec[4]
        name = "{}::{}::{}".format(contig, rec[1], ref) 
        annotations = rec[7]

        annotation_string = "ANN="

        #find snpEFF annotations

        try:
            #get only snpEFF anotations
            idx = annotations.index(annotation_string)
            snpeff_annotations = annotations[idx + 4 : ]
        except ValueError:
            snpeff_annotations = annotations

        annotation_regex = "[ATCGN]" + "\|([^|]*)"*4 + \
        "\|[^|]*"*4 + "\|([^|]*)\|([^|]*)\|"
 
        #the annotations field can contain multiple comma deliminated
        #entries

        annots = []

        for entry in snpeff_annotations.split(","):
            
            matches = re.search(annotation_regex, entry) 
  
            if matches:
                matched_entries = matches.groups()
                if annots:
                    annots = [",".join(i) for i in zip(annots, matched_entries)]
                else:
                    annots = matched_entries

        if not annots:
            annots = ["NA"] * 4 
        
        # get strand seperately
        strand_regex = "\|strand_([-+])\|"

        strands = re.search(strand_regex, snpeff_annotations)

        if strands:
            strand = strands.groups()[0]
        else:
            strand = "."

        if attrs is not None:
            # get use defined attributes
            attrs_annots = []
            for word in attrs:
                re_word = word + "=([^;]*);"
                re_matches = re.search(re_word, annotations)
                if re_matches:
                    rematched_entries = re_matches.groups()[0]
                    attrs_annots.append(rematched_entries)

            if not attrs_annots:
               attrs_annots = ["NA"] * attrs_len
       
            print(contig, start, end, ref, alt, name, \
                    strand, 
                    "\t".join(annots), \
                    annotations, \
                    "\t".join(attrs_annots), sep = "\t")
        else: 
            print(contig, start, end, ref, alt, name, \
                    strand, 
                    "\t".join(annots), \
                    annotations, \
                    sep = "\t")
        

def main():
    
    parser = argparse.ArgumentParser(description="""extract out records as
    bed Extended with annotation information """) 
    parser.add_argument('-v',
                          '--input_vcf',
                          help ='input vcf file)',
                       required = True),
    parser.add_argument('-a',
                          '--attrs',
                          nargs = '*',
                          help ='list of additional regexes to parse from vcf)',
                       required = False)
    args=parser.parse_args()
    
    gzopen = lambda f: gzip.open(f, 'rt') if f.endswith(".gz") else open(f)
    
    vcf = gzopen(args.input_vcf)
    attrs = args.attrs 
    
    vcf_to_bed(vcf, attrs)

    vcf.close()

if __name__ == '__main__': main()




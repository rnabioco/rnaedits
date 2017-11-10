#!/usr/bin/env python
import re, sys
import argparse

def gtf_to_bed(input, attributes, cols, output):
    
    for index, line in enumerate(input):
            
        if line.startswith("#"):
            continue
        
        attribute_list = []
        
        for attr in attributes:
            pat = '\\b' + attr + ' \"([^;.]+)\";'
            m = re.search(pat, line.split("\t")[8])
            if m is None:
                attribute_list.append("NA")
                continue
            
            else:
                attribute_list.append(m.group(1))
        
        if cols == "NA":
        
            pass
        
        else:
            
            for i in cols:
                
                try:
                    attribute_list.append(line.split("\t")[i - 1].rstrip('\n'))
                
                except IndexError:
                    sys.exit("Error column requested does not exist at line %d" % (index + 1))
                
        chrom = line.split("\t")[0]
        start = int(line.split("\t")[3]) - 1
        stop = int(line.split("\t")[4])
        strand = line.split("\t")[6]
        
        out_attrs = '\t'.join(str(e) for e in attribute_list)
        
        output.write('%s\t%d\t%d\t.\t.\t%s\t%s\n' % 
                (chrom, start, stop, strand, out_attrs)) 
                

if __name__ == '__main__':
    parser=argparse.ArgumentParser(description="""convert gtf to bed, while keeping attributes of interest""")
    parser.add_argument('-i','--input_gtf', help="""input gtf file""", required = True)
    parser.add_argument('-a','--attributes' , help="""list of attributes to keep (appended to the 7th to nth columns)""",required = True, nargs='+')
    parser.add_argument('-c','--cols' , help="""additional columns to keep (tab) deliminated""",required = False, nargs='+', type=int)
    parser.add_argument('-o','--output_bed' , help="""output bed with correct coordinates""",required = True)
    
    args = parser.parse_args()
    
    input_gtf = open(args.input_gtf, 'r')
    attributes = args.attributes
    output_bed = open(args.output_bed, 'w')
    cols=args.cols
    
    if cols:
        gtf_to_bed(input_gtf, attributes, cols, output_bed)
    else:
        gtf_to_bed(input_gtf, attributes,"NA", output_bed)
    
    output_bed.close() 

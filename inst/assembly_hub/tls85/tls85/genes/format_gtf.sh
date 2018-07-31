#! /usr/bin/env bash
# https://www.biostars.org/p/135790/

gtf="fixed.all95.gtf"
genome="chrom_sizes.txt"

gtfToGenePred -infoOut=infoOut.txt -genePredExt $gtf ${gtf/.gtf/.gp}

# Check the genePred output is valid
genePredCheck ${gtf/.gtf/.gp}

# Convert genePred format to BED format
genePredToBed ${gtf/.gtf/.gp} stdout | sort -k1,1 -k2,2n > ${gtf/.gtf/.bed}

# Convert BED to bigBed
# extraIndex required for position/search
bedToBigBed -type=bed12 -extraIndex=name ${gtf/.gtf/.bed} $genome ${gtf/.gtf/.bb}

# Required for indexing step
grep -v "^#" infoOut.txt | \
    awk '{printf "%s\t%s,%s,%s,%s,%s\n",$1,$2,$3,$8,$9,$10}' > ${gtf/.gtf/.nameIndex.txt}

# Create index for position/search function in browser
ixIxx ${gtf/.gtf/.nameIndex.txt} ${gtf/.gtf/.nameIndex.ix} ${gtf/.gtf/.nameIndex.ixx}

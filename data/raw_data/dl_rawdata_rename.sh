#! /usr/bin/env


# make array with sra id numbers from metadata

# download data from SRA database
for fq in fastqs
do echo $fq
  fastq-dump --gzip $fq
done

# rename files and move to requiste directory

cut -f 1,2,3 ../../docs/geo_formatted_metadata.txt \
    | awk -F'\t' 'system("mv " $1 " " $2/$3".fastq.gz")' 



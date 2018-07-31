#! /usr/bin/env bash
# http://genomewiki.ucsc.edu/index.php/RepeatMasker#Composite_RMSK_trackDb.txt
input_file=$1

mkdir rmskClass
sort -k12,12 $input_file | splitFileByColumn -ending=tab -col=12 -tab stdin rmskClass

for T in SINE LINE LTR DNA Simple Low_complexity Satellite Simple_repeat
do
 ./toBed6+10.pl rmskClass/${T}*.tab | sort -k1,1 -k2,2n > rmsk.${T}.bed
done
./toBed6+10.pl rmskClass/rRNA.tab | sort -k1,1 -k2,2n > rmsk.RNA.bed

for T in DNA LINE LTR Low_complexity RNA SINE Satellite Simple Simple_repeat
do
  bedToBigBed -tab -type=bed6+10 -as=rmsk16.as rmsk.${T}.bed \
   ../genes/chrom_sizes.txt rmsk.${T}.bb
done

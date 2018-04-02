#!/usr/bin/env bash

##################
# Copy assembly hub bigWigs files into organized directories
##################


#################
# Animal phenotype data
# Animal_number State
# 53  Ar
# 54  Ar
# 57  Ar
# 58  Ar
# 60  Ar
# 61  Ar
# 164 Ent
# 55  Ent
# 59  Ent
# 62  Ent
# 64  Ent
# 66  Ent
# 153 IBA
# 163 IBA
# 69  IBA
# 70  IBA
# 78  IBA
# 65  LT
# 67  LT
# 68  LT
# 72  LT
# 74  LT
# 29  SA
# 30  SA
# 31  SA
# 32  SA
# 33  SA
# 34  SA
# 77  SpD
# 79  SpD
# 80  SpD
# 83  SpD
# 84  SpD
#################


Ar=(53 54 57 58 60 61)
Ent=(164 55 59 62 64 66)
IBA=(153 163 69 70 78)
LT=(65 67 68 72 74)
SA=(29 30 31 32 33 34)
SpD=(77 79 80 83 84)

###############
# copy bams
###############

og_bg_dir="../../../../../data/bigwigs/"


for file in "${og_bg_dir}"*/*.bw
do echo "copying file "$file
    cp -u $file ./
done

states=("SA" "Ar" "Ent" "LT" "IBA" "SpD")

for state in ${states[@]}
do echo "making directory for state "$state
    mkdir -p $state
done

for file in ${SA[@]}
  do echo "moving SA files " *[MBH]$file*
      files=$(echo *[MBH]$file*);
      for good_file in ${files[@]}
        do mv $good_file SA/
        done
    done


for file in ${Ar[@]}
  do echo "moving Ar files " *[MBH]$file*
      files=$(echo *[MBH]$file*);
      for good_file in ${files[@]}
        do mv $good_file Ar/
        done
    done

for file in ${Ent[@]}
  do echo "moving Ent files " *[MBH]$file*
      files=$(echo *[MBH]$file*);
      for good_file in ${files[@]}
        do mv $good_file Ent/
        done
    done

for file in ${LT[@]}
  do echo "moving LT files " *[MBH]$file*
      files=$(echo *[MBH]$file*);
      for good_file in ${files[@]}
        do mv $good_file LT/
        done
    done

for file in ${IBA[@]}
  do echo "moving IBA files " *[MBH]$file*
      files=$(echo *[MBH]$file*);
      for good_file in ${files[@]}
        do mv $good_file IBA/
        done
    done

for file in ${SpD[@]}
  do echo "moving SpD files " *[MBH]$file*
      files=$(echo *[MBH]$file*);
      for good_file in ${files[@]}
        do mv $good_file SpD/
        done
    done


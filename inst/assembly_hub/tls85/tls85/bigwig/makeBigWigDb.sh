#!/usr/bin/env bash

# set directory to directory where bash script is located
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $DIR

#########
# general track settings
#########


echo "
track BigWigs
superTrack on
group martin_data
shortLabel bigWigs
longLabel bigWigs from RNA-Seq of animals at varying stages of hibernation

    track rnaseq_bigwigs
    parent BigWigs
    compositeTrack on
    shortLabel RNA-Seq
    longLabel bigWigs from RNA-Seq alignment data
    type bigWig
    autoScale on
    windowingFunction mean
    visibility hide 
    maxHeightPixels 100:30:8
    subGroup1 state Hibernation_State A_SA=Summer_Active B_IBA=Interbout_Arousal C_Ent=Entrance D_LT=Late_Torpor E_Ar=Arousal F_SpD=Spring_D
    subGroup2 region Brain_Region M=Medulla B=Brain_the_rest H=Hypothamalus
    subGroup3 rep Replicate R1=rep1 R2=rep2 R3=rep3 R4=rep4 R5=rep5
    subGroup4 strand Strand pos=Positive neg=Negative
    dimensions dimX=region dimY=state dimA=strand
    dimensionXchecked Medulla
    sortOrder strand=- state=+
"

########
# per sample bigwig settings
########

states=("SA" "IBA" "Ent" "LT" "Ar" "SpD")
colors=("255,0,0" "67,205,128" "155,48,255" "25,25,112" "0,0,255" "255,165,0")
prefix_for_state=("A_" "B_" "C_" "D_" "E_" "F_")
#set counter to grab color by index
color_idx=0

for state in ${states[@]}
do 
   
    #set counters to follow replicates
    which_rep=0
    med_rep=1
    hypo_rep=1
    fore_rep=1
    #get color value
    color=${colors[$color_idx]}
    
    #get prefix
    prefix=${prefix_for_state[$color_idx]}
    
    color_idx=$((color_idx + 1))

    for file in $state/*.bw
    do 
        med=$(echo $file | grep -c "M[0-9]")
        hypo=$(echo $file | grep -c "H[0-9]")
        fore=$(echo $file | grep -c "B[0-9]")
        pos_neg=$(echo $file | grep -c "pos")

        if [ $med -eq 1 ]
        then
            region="M"
            long_region="Medulla"
            which_rep=$med_rep
            if [ $pos_neg -eq 1 ]
            then
                med_rep=$((med_rep + 1))
            fi

        elif [ $hypo -eq 1 ]
        then
            region="H"
            long_region="Hypothalamus"
            which_rep=$hypo_rep
            if [ $pos_neg -eq 1 ]
            then
                hypo_rep=$((hypo_rep + 1))
            fi
                
        elif [ $fore -eq 1 ]
        then    
            region="B"
            long_region="Brain_the_rest"
            which_rep=$fore_rep
            if [ $pos_neg -eq 1 ]
            then
                fore_rep=$((fore_rep + 1))
            fi
        
        else 
            echo "Filename not found to be medulla hypothamalus or brain the rest"  >&2
            exit 1
        
        fi


        if [ $pos_neg -eq 1 ]
        then 
            strand="pos"
        else
            strand="neg"
        fi
        
        trackname=$(basename $file .bw) 

        echo "track $trackname"
        echo "parent rnaseq_bigwigs"
        echo "shortLabel $trackname"
        echo "bigDataUrl http://amc-sandbox.ucdenver.edu/User33/Martin/tls85/tls85/bigwig/$file"
        echo "longLabel  $trackname  from region $long_region in state $state replicate $which_rep"
        echo "subGroups state=$prefix$state region=$region replicate=R$which_rep" strand=$strand
        echo "maxHeightPixels 30:30:10"
        echo "color $color"
        echo "type bigWig"
        echo ""
    
   done
done

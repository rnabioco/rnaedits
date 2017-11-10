#! /usr/bin/env bash

# download required reference databases

# make organized directories
mkdir -p ensembl85 denovo_annotation

ens="ensembl85"
denovo="denovo_annotation"

# Ensembl gtf
wget -P $ens ftp://ftp.ensembl.org/pub/release-85/gtf/ictidomys_tridecemlineatus/Ictidomys_tridecemlineatus.spetri2.85.gtf.gz

# denovo transcriptome
wget -P $denovo http://amc-sandbox.ucdenver.edu/User33/Martin/fixed.all95.gtf

# genome fasta
wget -P $ens ftp://ftp.ensembl.org/pub/release-85/fasta/ictidomys_tridecemlineatus/dna/Ictidomys_tridecemlineatus.spetri2.dna.toplevel.fa.gz

# index fasta
gunzip $ens"/Ictidomys_tridecemlineatus.spetri2.dna.toplevel.fa.gz"
samtools faidx $ens"/Ictidomys_tridecemlineatus.spetri2.dna.toplevel.fa"

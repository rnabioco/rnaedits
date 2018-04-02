#! /usr/bin/env bash

gatk="/vol3/home/riemondy/Projects/Martin/data/vcf/variant_allele_counts_by_strand/A_G_alleles"


tabix -p vcf $gatk"/filtered_select_variants.vcf.gz"

cp $gatk"/filtered_select_variants.vcf.gz" \
    $gatk"/filtered_select_variants.vcf.gz.tbi" \
    ./



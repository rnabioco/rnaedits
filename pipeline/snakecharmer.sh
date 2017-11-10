#!/usr/bin/env bash
#BSUB -J RNAedits 
#BSUB -o logs/snakemake_%J.out
#BSUB -e logs/snakemake_%J.err
#BSUB -R "select[mem>4] rusage[mem=4] " 
#BSUB -q normal

set -o nounset -o pipefail -o errexit -x

args=' -q rna -o {log}.out -e {log}.err -J {params.job_name} -R "
{params.memory} span[hosts=1] " -n {threads} '

snakemake --drmaa "$args" \
    --snakefile Snakefile \
    --jobs 75 \
    --resources all_threads=75 \
    --latency-wait 50 \
    --rerun-incomplete \
    --configfile config_denovo.yaml 
  

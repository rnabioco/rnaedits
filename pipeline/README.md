## Snakemake pipeline

This snakemake pipeline was used to identify rna editing sites and
was executed on a cluster using an LSF DRMAA scheduler. This cluster
uses modules to manage software and to run this pipeline calls
to modules will need to be removed, and the the relevant software placed
into PATH for proper execution.

This pipeline has the following dependencies: 

```bash
snakemake 

GATK/3.8
picard/2.7.0
samtools/1.5
htslib/1.5
bcftools/1.5
cutadapt/1.16
fastqc/0.11.5
bwa/0.7.16a-r1181
subread/1.6.2
vcflib/7e3d806  
star/2.5.2a
bedtools/2.26.0

# to run python scripts
python
python3
pandas
pysam

# to run R scripts
R
dplyr
purrr
tidyr
readr
stringr
ggplot2
RColorBrewer
edgeR
ComplexHeatmap
VennDiagram
gridExtra
cowplot
```

## Download annotations

```bash
cd dbases
bash dl_data.sh
```

## Download raw data

The raw data contains 90 paired end libraries. 
```bash
cd data/raw_data
bash dl_data.sh
```

## Build C++ scripts

These scripts are used by the hyperediting pipeline
```bash
cd src
make
```

## Edit config file


The config.yaml file contains hardcoded paths to various tools and 
annotations that will need to modified to match your respective paths. 

## Test pipeline

To show output of main pipeline (~2000 jobs)
```bash
cd pipeline
snakemake -npr
```

To show output of hyperediting pipeline (~22000 jobs)
```bash
snakemake -npr -s hyperediting.snake
```


## Run pipeline
To execute both pipelines (~ 14 days on a cluster using up to 70 cores)

```
bsub < snakecharmer.sh
```


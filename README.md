## rnaedits

Pipeline and scripts for detecting A-to-I RNA editing used in the publication:

> Riemondy, K. A., Gillen, A. E., White, E. A., Bogren, L. K., Hesselberth, J. R., & Martin, S. L. (2018). 
*Dynamic temperature-sensitive A-to-I RNA editing in the brain of a heterothermic mammal during hibernation.* **RNA**. [https://doi.org/10.1261/rna.066522.118]()

## UCSC Trackhub

BigWigs and editing sites are viewable as an assembly hub in the UCSC
genome browser. 

To load assembly hub:

1) go to [My Data](https://genome.ucsc.edu/cgi-bin/hgHubConnect) on the
UCSC website

2) paste the following link into the `My Hubs` URL input box. 

`http://amc-sandbox.ucdenver.edu/User33/Martin/tls85/hub.txt`  

3) click add hub

## Pipelines used to detect RNA editing

A snakemake pipeline was used to run a pipeline to detect A-to-I editing
events. This pipeline requies a cluster to run all of the steps, and is
provided in the `pipeline` directory. 

## Rmarkdown documents used to generate figures

Rmarkdown documents are located in the `results` directory. A
snakemake pipeline is provided to run the rmarkdown (`results/notebook`)



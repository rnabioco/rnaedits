## Hyperediting Detection
library(dplyr)
library(tidyverse)
library(stringr)
library(tidyr)
library(readr)

args<-commandArgs(TRUE)

data_dir <- args[1]

alleles  <- list(
              "A_to_C",
              "A_to_G",
              "A_to_T",
              "C_to_A",
              "C_to_G",
              "C_to_T",
              "G_to_A",
              "G_to_C",
              "G_to_T",
              "T_to_A",
              "T_to_C",
              "T_to_G")

# read and parse hyperedited reads/sites
read_hyperedits <- function(allele, type = "reads"){
  # first load in bed files of individual hyperedited reads
  files <- dir(file.path(data_dir, "hyperedits", "bed",
                     type, "filtered", allele),
           pattern = "*.bed", 
           recursive = T,
           full.names = T)
  
  dat <- map(files, ~read_tsv(.x, col_names = c("chrom",
                                                "start",
                                                "end",
                                                "id",
                                                "mismatches",
                                                "strand"), 
                              col_types = c("cddccc")))
  file_ids <- basename(files)
  regions <- str_extract(files, "medulla|brainrest|hypothalamus")
  file_names <- str_c(file_ids, ":::", regions)
  names(dat) <- file_names
  
  dat <- bind_rows(dat, .id = "name")
  # exit if no sites/reads exist
  if (nrow(dat) == 0){ return(NULL) }
  
  dat <- separate(dat, name, c("fileid", "region"), ":::")
  dat <- mutate(dat,
                animal_number = str_split(fileid, "_", simplify = T)[, 1],
                animal_number = str_replace(animal_number, "^[0-9]+-", ""))
  #annotate with State
  
  pdata <- readr::read_tsv(file.path(docs_dir, "BrainRegionRNAseqLibMetadata.txt"))
  pdata <- gather(pdata, region, sample, Forebrain, Hypothalamus, Medulla)
  pdata <- mutate(pdata, 
                  animal_number = str_split(sample, "_", simplify = T)[, 1])
  simple_pdata <- pdata %>% 
    dplyr::select(State, animal_number)
  
  dat <- left_join(dat, 
                   simple_pdata, 
                   by ="animal_number")
  
  dat <- mutate(dat, 
                   region = ifelse(region == "brainrest",
                                   "Forebrain",
                                   ifelse(region == "medulla",
                                          "Medulla",
                                          ifelse(region == "hypothalamus",
                                                 "Hypothalamus",
                                                 NA))))
  
  dat
}

edit_dat <- map(alleles, 
                ~read_hyperedits(.x, type = "edits"))

names(edit_dat) <- alleles
edit_dat <- bind_rows(edit_dat, .id ="allele")

edit_dat %>% 
  dplyr::group_by(allele, 
                  strand, 
                  chrom, 
                  start, 
                  end) %>% 
  mutate(supporting_reads = length(unique(id))) %>% 
  ungroup() -> edit_dat

dplyr::filter(edit_dat, 
              supporting_reads > 1) -> edit_dat_1  

hyper_edits <- dplyr::filter(edit_dat_1, allele == "A_to_G") %>% 
  mutate(site_id = paste(chrom, 
                         start, 
                         ifelse(strand == "+", 
                                "A",
                                "T"), 
                         sep = "::"))


#vcf format
#chrom  pos (1bases)  ID  REF Alt QUAL  FILTER  INFO
# just make a pseudo-vcf, only need vcf format for get_pileup.py script
vcf_out <- dplyr::select(hyper_edits,
                         chrom, 
                         end, #1bases
                         site_id,
                         strand) %>% 
  unique()


bed_out <-  dplyr::mutate(hyper_edits,
                          name = ifelse(strand == "+",
                                        "strand_+",
                                        "strand_-")) %>% 
  dplyr::select(chrom, start, end, #1bases
                name, site_id, strand ) %>% 
  dplyr::mutate_if(is.double, as.integer) %>% 
  ungroup() %>% 
  unique()

vcf_out <- dplyr::mutate(vcf_out,
                         `#CHROM` = chrom,
                         POS = as.integer(end),
                         ID = ".",
                         REF = ifelse(strand == "+",
                                      "A",
                                      "T"),
                         ALT = ifelse(strand == "+",
                                      "G",
                                      "C"),
                         QUAL = 0.0,
                         FILTER = "PASS",
                         INFO = "AC=1",
                         FORMAT =  "GT:AD:DP:GQ:PL",
                         Dummy_sample = "1/1:0,2:2:6:1,6,0") %>% 
  ungroup() %>% 
  dplyr::select(-c(strand, chrom, end, site_id))

write_tsv(bed_out, file.path(data_dir, "hyperedits", "hyperedited_sites.bed"), col_names = T)
write_tsv(vcf_out, file.path(data_dir, "hyperedits", "hyperedited_sites.vcf"), col_names = T)





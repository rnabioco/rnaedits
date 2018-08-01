# libraries
library(dplyr)
library(purrr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(edgeR)
library(ComplexHeatmap)
library(VennDiagram)
library(gridExtra)
library(cowplot)
library(valr)
library(flexclust) # for kmeans++ clustering
library(kentr) # to get sequences using htslib
library(topGO)
library(scales)

# colors
state_cols <-  c(
  SA = rgb(255, 0, 0, maxColorValue=255),
  IBA = rgb(67, 205, 128, maxColorValue=255),
  Ent = rgb(155, 48, 255, maxColorValue=255),
  LT = rgb(25, 25, 112, maxColorValue=255),
  Ar = rgb(0, 0, 255, maxColorValue=255),
  SpD = rgb(255, 165, 0, maxColorValue=255)
)

state_order = c(
  "SA", "IBA", "Ent", "LT", "Ar", "SpD"
)

#colors for region
region_cols <- c(
  Medulla = "#4DAF4A",
  Hypothalamus = "#377EB8",
  Forebrain =  "#E41A1C",
  Cerebrum = "#E41A1C"
)

# region order
region_order <- c("Medulla", "Hypothalamus", "Forebrain")


color_fxn <- colorRampPalette(brewer.pal(9, "Spectral"))

#### Paths ####

project_dir <- path.expand("~/Projects/publication_repos/rnaedits/")
data_dir <- file.path(project_dir, "data")
results_dir <- file.path(project_dir, "results")
docs_dir <- file.path(project_dir, "docs")
db_dir <- file.path(project_dir, "dbases")

# vector of figure paths
main_figs <- file.path(results_dir, "Figures", paste0("Figure_", 1:6))
sup_figs <- file.path(results_dir, "Figures", paste0("Sup_Figure_", 1:17))
all_figs <- c(main_figs, sup_figs)

for(i in seq_along(all_figs)){
  if(!dir.exists(all_figs[i])){
    dir.create(all_figs[i], showWarnings = FALSE, recursive = TRUE)
  } 
}

figs_dir <-  file.path(results_dir, "Figures") %>%
  dir(pattern = "^Figure_[1-6]$",
      include.dirs = TRUE,
      full.names = T)

sfigs_dir <-  file.path(results_dir, "Figures") %>%
  dir(pattern = "Sup_Figure_[0-9]$",
      include.dirs = TRUE,
      full.names = T)

sfigs_dir <- c(sfigs_dir,
               file.path(results_dir, "Figures") %>%
                 dir(pattern = "Sup_Figure_[1][0-9]$",
                     include.dirs = TRUE,
                     full.names = T)
)
##### Functions ####

#' When writing out excel workbooks using openxlsx::write.xlsx()
#' this function will set the class attributes for a column, which
#' enforces a column type in the resulting xlsx file. 
#' Useful for avoid gene names being clobbered to dates and 
#' setting scientific number formatting

set_xlsx_class <- function(df, col, xlsx_class){
  for(i in seq_along(col)){
    class(df[[col[i]]]) <- xlsx_class
  }
  df
}



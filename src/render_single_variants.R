args<-commandArgs(TRUE)


render_variant <- function(allele, indir, outdir, RMD_doc, variant_type){
  allele_path <- paste0(indir, allele, "/counts/")
  variant_info_path <- paste0(indir, allele, "/filtered_select_variants.bed.gz")

  ref_base <- stringr::str_split(allele, "_", simplify = T)[1]
  alt_base <- stringr::str_split(allele, "_", simplify = T)[2]
  message(allele_path)
  message(variant_info_path)
  message(paste("working on ", allele, " to ", ref_base, " ", alt_base))
  rmarkdown::render(RMD_doc,
                    output_file = paste0(outdir, allele, "/", allele, ".html"),
                    output_dir = paste0(outdir, allele),
                    clean = FALSE,
                    params = list(
                      allele_counts = allele_path,
                      variant_info = variant_info_path,
                      ref = ref_base,
                      alt = alt_base,
                      prefix = file.path("..", variant_type),
                      dirname = allele
                    ))
  }

Sys.setenv(RSTUDIO_PANDOC="/beevol/home/riemondy/bin/pandoc/")
render_variant(args[1], args[2], args[3], args[4],  args[5])

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(GGally)
  library(rmarkdown)
})

print(abundances)
print(rmdtemplate)
print(outhtml)

abundances <- readRDS(abundances)

plot_abundance_report <- function(abundances, output_file, output_dir = "./",
                                  output_format = "html_document", 
                                  rmd_template = NULL,
                                  knitr_progress = FALSE, ...){
  output_report <- file.path(output_dir, basename(output_file))
  output_rmd <- file.path(output_dir,
                          paste0(tools::file_path_sans_ext(basename(output_file)),
                                 ".Rmd"))
  template_file <- rmd_template
  file.copy(from = template_file, to = output_rmd, overwrite = TRUE)
  
  args <- list(...)
  args$input <- output_rmd
  args$output_format <- output_format
  args$output_file <- output_file
  args$quiet <- !knitr_progress

  output_file <- do.call("render", args = args)
  invisible(output_file)
}

plot_abundance_report(abundances = abundances, 
                      output_file = basename(outhtml), 
                      output_dir = dirname(outhtml),
                      rmd_template = rmdtemplate)

sessionInfo()
date()

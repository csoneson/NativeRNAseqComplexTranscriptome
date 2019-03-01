args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(gff)
print(gtf)

suppressPackageStartupMessages({
  library(rtracklayer)
  library(readr)
  library(withr)
})

## Filter out exons that can not be handled by gffcompare
x <- readr::read_tsv(gff, col_names = FALSE, col_types = "cccddcccc")
dim(x)
message("Excluding the following lines:")
x[x$X5 - x$X4 >= 30000, ]
x <- x[x$X5 - x$X4 < 30000, ]
dim(x)
withr::with_options(c(scipen = 100), 
                    write.table(x, file = gsub("\\.gff$", ".fixed.gff", gff), 
                                quote = FALSE, sep = "\t", row.names = FALSE, 
                                col.names = FALSE))

x <- rtracklayer::import(gsub("\\.gff$", ".fixed.gff", gff))
x$transcript_id <- as.character(x$group)
x$group <- NULL
rtracklayer::export(x, gtf)

date()
sessionInfo()
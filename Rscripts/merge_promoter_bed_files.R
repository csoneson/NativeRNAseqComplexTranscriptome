args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(GenomeInfoDb)
})

bedsin <- strsplit(bedsin, ",")[[1]]

print(bedsin)
print(bedout)

## Read all input bed files and merge into one GRanges object
suppressWarnings(bedall <- do.call(c, lapply(bedsin, function(b) rtracklayer::import(b))))

## Reduce the ranges in the merged GRanges object
bedrr <- GenomicRanges::reduce(bedall)

## Change the seqlevels style to Ensembl
GenomeInfoDb::seqlevelsStyle(bedrr) <- "Ensembl"

## Write reduced promoter regions to bed file
rtracklayer::export(bedrr, con = bedout, format = "bed")

date()
sessionInfo()
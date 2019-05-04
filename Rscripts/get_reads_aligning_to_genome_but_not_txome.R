args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(nbrreadsrds)
print(outdir)

## Extract read names for reads mapping only to genome (not transcriptome)

rd <- readRDS(nbrreadsrds)
rdg <- rd$genomebams
rdt_p0.99 <- rd$txomebams_p0.99

stopifnot(names(rdg) %in% names(rdt_p0.99),
          names(rdt_p0.99) %in% names(rdg))

for (nm in names(rdg)) {
  readnames <- setdiff(rdg[[nm]]$allAlignments$read,
                       rdt_p0.99[[nm]]$allAlignments$read)
  write.table(readnames, 
              file = paste0(outdir, "/", nm, 
                            "_reads_aligning_to_genome_but_not_txome.txt"),
              row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
}

saveRDS(NULL, file = paste0(outdir, "/reads_aligning_to_genome_but_not_txome.rds"))
date()
sessionInfo()

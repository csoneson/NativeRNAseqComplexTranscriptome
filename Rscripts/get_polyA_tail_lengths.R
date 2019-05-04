args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

## Estimate the length of the polyA tail

print(fast5dir)
print(ncores)
print(outdir)
print(csvfilename)

suppressPackageStartupMessages({
  library(tailfindr)
})

df <- find_tails(fast5_dir = fast5dir,
                 save_dir = outdir,
                 csv_filename = csvfilename,
                 num_cores = ncores)

saveRDS(df, file = paste0(outdir, "/", csvfilename, ".rds"))
date()
sessionInfo()


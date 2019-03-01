args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

## Plot number of reads per equivalence class, as function of number of
## transcripts per equivalence class

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

print(salmondir)
print(outrds)

salmondirs <- list.files(salmondir, full.names = TRUE)
salmonfiles <- paste0(salmondirs, "/aux_info/eq_classes.txt")
names(salmonfiles) <- basename(salmondirs)
salmonfiles <- salmonfiles[file.exists(salmonfiles)]
salmonfiles

eqcl <- lapply(salmonfiles, function(f) {
  x <- readLines(f)
  n_tr <- as.numeric(x[1])  ## Total number of transcripts
  n_eq <- as.numeric(x[2])  ## Total number of equivalence classes
  tx_id <- x[3:(n_tr + 2)]  ## Transcript IDs
  quants <- x[(n_tr + 3):length(x)]  ## Characteristics of equivalence classes
  
  ## Split equivalence class characteristics. Each element of the list corresponds
  ## to one equivalence class, and lists its number of transcripts, the
  ## transcripts IDs and the total number of reads
  do.call(rbind, lapply(quants, function(w) {
    tmp = strsplit(w, "\\\t")[[1]]
    nbr_tx = as.numeric(tmp[1])
    data.frame(nbr_tx = nbr_tx,
               tx_ids = paste0(tx_id[as.numeric(tmp[2:(1 + nbr_tx)]) + 1], collapse = ","),
               count = as.numeric(tmp[length(tmp)]),
               stringsAsFactors = FALSE)
  }))
})

for (nm in names(eqcl)) {
  eqcl[[nm]]$sample <- nm
}

eqcl <- do.call(rbind, eqcl)

pdf(gsub("rds$", "pdf", outrds))
ggplot(eqcl, aes(x = nbr_tx, y = count)) + geom_point(size = 0.5, alpha = 0.25) + 
  facet_wrap(~ sample) + theme_bw() + scale_y_log10() + 
  xlab("Number of transcripts in equivalence class") + 
  ylab("Number of reads assigned to equivalence class")
dev.off()

saveRDS(NULL, file = outrds)

sessionInfo()
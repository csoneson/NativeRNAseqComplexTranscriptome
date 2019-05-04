args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

sjdbfiles <- strsplit(sjdbfiles, ",")[[1]]

print(sjdbfiles)
print(min_unique)
print(outtxt)

suppressPackageStartupMessages({
  library(dplyr)
})

x <- do.call(dplyr::bind_rows, lapply(sjdbfiles, function(f) {
  read.delim(f, header = FALSE, as.is = TRUE) %>%
    setNames(c("chromosome", "start", "end", "strand", "motif", 
               "annotated", "nbr_unique", "nbr_mm", "max_overhang"))
})) %>%
  dplyr::group_by(chromosome, start, end, strand) %>%
  dplyr::summarize(motif = motif[1], 
                   annotated = annotated[1],
                   nbr_unique = sum(nbr_unique),
                   nbr_mm = sum(nbr_mm),
                   max_overhang = max(max_overhang)) %>%
  dplyr::filter(nbr_unique >= min_unique) %>%
  dplyr::select(chromosome, start, end, strand, motif, 
                annotated, nbr_unique, nbr_mm, max_overhang)

write.table(x, file = outtxt, col.names = FALSE, row.names = FALSE,
            quote = FALSE, sep = "\t")

date()
sessionInfo()

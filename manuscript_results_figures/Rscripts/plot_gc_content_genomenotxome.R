args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(stringr)
  library(Biostrings)
  library(data.table)
})
source("manuscript_results_figures/Rscripts/remap_sample_names.R")

datasets <- strsplit(datasets, ",")[[1]]
names(datasets) <- datasets
conditions <- strsplit(conditions, ",")[[1]]

print(datasets)
print(conditions)
print(outrds)

fastqs <- do.call(dplyr::bind_rows, lapply(datasets, function(ds) {
  fastqdir <- paste0(ds, "/FASTQ", 
                     ifelse(ds %in% c("RNA001", "HEK293RNA"), "dna", ""))
  fastqfiles <- list.files(fastqdir, pattern = "(FASTQ|fastq)\\.gz$", 
                           full.names = TRUE)
  names(fastqfiles) <- gsub("\\.(FASTQ|fastq).gz", "", basename(fastqfiles))
  fastqfiles <- fastqfiles[!grepl("_orig", names(fastqfiles))]
  fastqfiles <- fastqfiles[names(fastqfiles) %in% 
                             sample_annotation$sample_orig[sample_annotation$condition %in% conditions]]
  fastqfiles
  
  genomenotxomedir <- paste0(ds, "/FASTQgenomenotxome")
  textfiles <- list.files(genomenotxomedir, pattern = "_reads_aligning_to_genome_but_not_txome.txt",
                          full.names = TRUE)
  names(textfiles) <- gsub("_reads_aligning_to_genome_but_not_txome.txt", "", basename(textfiles))
  textfiles <- textfiles[!grepl("_orig", names(textfiles))]
  textfiles <- textfiles[names(textfiles) %in% 
                           sample_annotation$sample_orig[sample_annotation$condition %in% conditions]]
  textfiles
  
  stopifnot(names(fastqfiles) %in% names(textfiles),
            names(textfiles) %in% names(fastqfiles))
  textfiles <- textfiles[match(names(fastqfiles), names(textfiles))]
  stopifnot(names(fastqfiles) == names(textfiles))
  
  do.call(dplyr::bind_rows, lapply(names(fastqfiles), function(nm) {
    f <- fastqfiles[nm]
    fastq <- fread(paste0("zcat ", f), sep = "\n", header = FALSE)$V1
    seq_idxs <- seq(2, length(fastq), by = 4)
    reads <- Biostrings::DNAStringSet(x = fastq[seq_idxs])
    names(reads) <- sapply(strsplit(fastq[seq_idxs - 1], " "), .subset, 1)
    gc <- letterFrequency(reads, "GC", as.prob = TRUE)
    gnt <- read.delim(textfiles[nm], header = FALSE, as.is = TRUE)$V1
    as.data.frame(gc) %>% dplyr::mutate(read = names(reads), 
                                        sample = remap[nm],
                                        dataset = remapds[ds]) %>%
      dplyr::mutate(genomenotxome = gsub("^@", "", read) %in% gnt)
  }))
}))

fastqs <- fastqs %>% 
  dplyr::mutate(readcat = c("Other reads", "Reads mapping to genome but not transcriptome")[genomenotxome + 1])
print(table(fastqs$sample, fastqs$readcat))

png(gsub("rds$", "png", outrds), width = 6, height = 3, unit = "in", res = 400)
ggplot(fastqs, aes(x = `G|C`, group = interaction(sample, readcat), color = readcat)) + 
  geom_line(stat = "density", size = 1.5) + 
  theme_bw() + 
  scale_color_manual(values = c("red", "blue"), name = "") + 
  xlab("GC content") + ylab("Density") + 
  theme(legend.position = "bottom")
dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()

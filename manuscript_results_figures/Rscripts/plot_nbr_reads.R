args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(plyr)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(stringr)
})
source("manuscript_results_figures/Rscripts/remap_sample_names.R")

datasets <- strsplit(datasets, ",")[[1]]
names(datasets) <- datasets
conditions <- strsplit(conditions, ",")[[1]]

print(datasets)
print(conditions)
print(outrds)

## ========================================================================== ##
## Total number of reads
## ========================================================================== ##
## Read information about the number of reads
reads <- lapply(datasets, function(ds) {
  read.delim(paste0(ds, "/output/", ds, "_nbr_reads.txt"),
             header = TRUE, as.is = TRUE) %>%
    dplyr::filter(!grepl("_orig", sample)) %>%
    dplyr::mutate(sample = remap[sample],
                  dataset = remapds[ds]) %>%
    dplyr::select(sample, dataset, nReads)
})

## Put everything together
reads <- do.call(dplyr::bind_rows, reads) %>% 
  dplyr::left_join(sample_annotation %>% dplyr::select(sample_remap, condition),
                   by = c("sample" = "sample_remap")) %>%
  dplyr::filter(condition %in% conditions)

## Write total number of reads to text file
write.table(reads, file = gsub("\\.rds$", "_nbrreads.txt", outrds),
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

## Plot total number of reads
ptot <- ggplot(
  reads %>% dplyr::mutate(sample = removeDatasetFromSample(sample, dataset)), 
  aes(x = sample, y = nReads/1e6), fill = "#777777") + 
  geom_bar(stat = "identity", position = "dodge") + 
  scale_y_continuous(expand = c(0, 0, 0.05, 0)) + 
  facet_grid(~ dataset, scales = "free_x", space = "free_x") + theme_bw() + 
  theme(legend.position = "none",
        strip.text = element_text(size = 8)) + 
  xlab("") + ylab("Total number of reads (Mio.)")

## ========================================================================== ##
## Read length
## ========================================================================== ##
## Read information about read lengths
readInfo <- do.call(dplyr::bind_rows, lapply(datasets, function(ds) {
  rd <- readRDS(paste0(ds, "/output/", ds, "_nbr_reads.rds"))
  rdf <- rd$fastqs
  rdf <- rdf[!grepl("_orig", names(rdf))]
  rdf <- rdf[names(rdf) %in% 
               sample_annotation$sample_orig[sample_annotation$condition %in% conditions]]
  do.call(dplyr::bind_rows, lapply(names(rdf), function(nm) {
    rdf[[nm]]$reads %>% dplyr::select(readLength, aveBaseQuality) %>%
      dplyr::mutate(sample = remap[nm], dataset = remapds[ds])
  }))
}))

## Write summary information about read lengths to text file
write.table(
  readInfo %>% dplyr::group_by(sample) %>%
    dplyr::summarize(mean_readlength = mean(readLength),
                     median_readlength = median(readLength),
                     min_readlength = min(readLength),
                     max_readlength = max(readLength)),
  file = gsub("\\.rds$", "_readlengths.txt", outrds),
  quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
)

## Plot read lengths
plength <- ggplot(readInfo %>% 
                    dplyr::mutate(replicate = sapply(strsplit(sample, "_"), 
                                                     .subset, 3)), 
                  aes(x = readLength, group = sample, alpha = as.factor(replicate))) + 
  geom_line(stat = "density", size = 1) + scale_x_log10() + 
  facet_wrap(~ dataset, nrow = 2, scales = "free_y") + theme_bw() + 
  xlab("Read length") + ylab("Density") + 
  scale_alpha_manual(values = c(1, 0.2, 0.68, 0.36, 0.84, 0.52)) + 
  theme(legend.position = "none")

## Plot mean quality distribution
pqual <- ggplot(readInfo %>% 
                  dplyr::mutate(replicate = sapply(strsplit(sample, "_"), 
                                                   .subset, 3)), 
                aes(x = aveBaseQuality, group = sample, alpha = as.factor(replicate))) + 
  geom_line(stat = "density", size = 1) + 
  facet_wrap(~ dataset, nrow = 2, scales = "free_y") + theme_bw() + 
  xlab("Average base quality per read") + ylab("Density") + 
  scale_alpha_manual(values = c(1, 0.2, 0.68, 0.36, 0.84, 0.52)) + 
  theme(legend.position = "none")

## Plot read length vs mean quality, for each sample
plengthqual <- 
  ggplot(readInfo, aes(x = readLength, y = aveBaseQuality)) + 
  geom_hex(bins = 100, aes(fill = stat(density))) + scale_x_log10() + 
  scale_fill_gradient(name = "", low = "bisque2", high = "darkblue") + 
  facet_wrap(~ sample) + theme_bw() + 
  xlab("Read length") + ylab("Average base quality")

## Make final plots
png(gsub("\\.rds$", "_nbrreads.png", outrds), width = 12, 
    height = 8, unit = "in", res = 400)
print(cowplot::plot_grid(
  ptot,
  cowplot::plot_grid(plength, pqual, nrow = 1, 
                     rel_widths = c(1, 1), labels = c("B", "C")),
  ncol = 1, rel_heights = c(1, 1), labels = c("A", ""))
)
dev.off()

png(gsub("\\.rds$", "_readlength_vs_quality.png", outrds), width = 11, 
    height = 11, unit = "in", res = 400)
print(plengthqual)
dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()

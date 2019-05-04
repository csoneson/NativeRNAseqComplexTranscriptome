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
  dplyr::filter(condition %in% conditions) %>% 
  dplyr::mutate(
    dataset = factor(dataset, levels = ds_order[ds_order %in% dataset])
  )

## Write total number of reads to text file
write.table(reads, file = gsub("\\.rds$", "_nbrreads.txt", outrds),
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

## Plot total number of reads
ptot <- ggplot(
  reads %>% dplyr::mutate(sample = removeDatasetFromSample(sample, dataset)), 
  aes(x = sample, y = nReads/1e6, fill = dataset)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  scale_y_continuous(expand = c(0, 0, 0.05, 0)) + 
  facet_grid(~ dataset, scales = "free_x", space = "free_x") + theme_bw() + 
  scale_fill_manual(values = ds_colors, name = "") + 
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
})) %>% 
  dplyr::mutate(
    dataset = factor(dataset, levels = ds_order[ds_order %in% dataset])
  )

## Write summary information about read lengths to text file
write.table(
  readInfo %>% dplyr::group_by(sample) %>%
    dplyr::summarize(mean_readlength = mean(readLength),
                     median_readlength = median(readLength),
                     min_readlength = min(readLength),
                     max_readlength = max(readLength),
                     nbr_reads_longer_than_15000 = length(which(readLength > 1.5e4))),
  file = gsub("\\.rds$", "_readlengths.txt", outrds),
  quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
)

## -------------------------------------------------------------------------- ##
## Plot read lengths
plength <- ggplot(readInfo %>% 
                    dplyr::mutate(replicate = sapply(strsplit(sample, "_"), 
                                                     .subset, 3)), 
                  aes(x = readLength, group = sample, 
                      color = dataset, alpha = as.factor(replicate))) + 
  geom_line(stat = "density", size = 1) + scale_x_log10() + 
  facet_wrap(~ dataset, nrow = 2, scales = "free_y") + theme_bw() + 
  xlab("Read length") + ylab("Density") + 
  scale_color_manual(values = ds_colors, name = "") + 
  scale_alpha_manual(values = c(1, 0.4, 0.76, 0.52, 0.88, 0.64)) + 
  theme(legend.position = "none")

## Plot read lengths, sqrt-transformed
plengthsqrt <- ggplot(readInfo %>% 
                        dplyr::mutate(replicate = sapply(strsplit(sample, "_"), 
                                                         .subset, 3)), 
                      aes(x = readLength, group = sample, 
                          color = dataset, alpha = as.factor(replicate))) + 
  geom_line(stat = "density", size = 1) + scale_x_sqrt(limits = c(0, 1.5e4)) + 
  facet_wrap(~ dataset, nrow = 2, scales = "free_y") + theme_bw() + 
  xlab("Read length") + ylab("Density") + 
  scale_color_manual(values = ds_colors, name = "") + 
  scale_alpha_manual(values = c(1, 0.4, 0.76, 0.52, 0.88, 0.64)) + 
  theme(legend.position = "none")

## Same plot but with x-axis on linear scale, and cut at 1.5e4
plengthlin <- ggplot(readInfo %>% 
                       dplyr::mutate(replicate = sapply(strsplit(sample, "_"), 
                                                        .subset, 3)), 
                     aes(x = readLength, group = sample, 
                         color = dataset, alpha = as.factor(replicate))) + 
  geom_line(stat = "density", size = 1) + 
  xlim(0, 1.5e4) + 
  facet_wrap(~ dataset, nrow = 2, scales = "free_y") + theme_bw() + 
  xlab("Read length") + ylab("Density") + 
  scale_color_manual(values = ds_colors, name = "") + 
  scale_alpha_manual(values = c(1, 0.4, 0.76, 0.52, 0.88, 0.64)) + 
  theme(legend.position = "none")

## -------------------------------------------------------------------------- ##
## Plot mean quality distribution
pqual <- ggplot(readInfo %>% 
                  dplyr::mutate(replicate = sapply(strsplit(sample, "_"), 
                                                   .subset, 3)), 
                aes(x = aveBaseQuality, group = sample, 
                    color = dataset, alpha = as.factor(replicate))) + 
  geom_line(stat = "density", size = 1) + 
  facet_wrap(~ dataset, nrow = 2, scales = "free_y") + theme_bw() + 
  xlab("Average base quality per read") + ylab("Density") + 
  scale_color_manual(values = ds_colors, name = "") + 
  scale_alpha_manual(values = c(1, 0.4, 0.76, 0.52, 0.88, 0.64)) + 
  theme(legend.position = "none")

## -------------------------------------------------------------------------- ##
## Plot read length vs mean quality, for each sample
plengthqual <- 
  ggplot(readInfo %>% dplyr::arrange(dataset, sample) %>%
           dplyr::mutate(sample = factor(sample, levels = unique(sample))), 
         aes(x = readLength, y = aveBaseQuality)) + 
  geom_hex(bins = 100, aes(fill = stat(density))) + scale_x_log10() + 
  scale_fill_gradient(name = "", low = "bisque2", high = "darkblue") + 
  facet_wrap(~ sample) + theme_bw() + 
  xlab("Read length") + ylab("Average base quality")

## Same plot but with x-axis on linear scale, and cut at 1.5e4
plengthquallin <- 
  ggplot(readInfo %>% dplyr::arrange(dataset, sample) %>%
           dplyr::mutate(sample = factor(sample, levels = unique(sample))), 
         aes(x = readLength, y = aveBaseQuality)) + 
  geom_hex(bins = 100, aes(fill = stat(density))) + 
  xlim(0, 1.5e4) + 
  scale_fill_gradient(name = "", low = "bisque2", high = "darkblue") + 
  facet_wrap(~ sample) + theme_bw() + 
  xlab("Read length") + ylab("Average base quality")

## Plot read length vs mean quality, for each sample, sqrt-transformed
plengthqualsqrt <- 
  ggplot(readInfo %>% dplyr::arrange(dataset, sample) %>%
           dplyr::mutate(sample = factor(sample, levels = unique(sample))), 
         aes(x = readLength, y = aveBaseQuality)) + 
  geom_hex(bins = 100, aes(fill = stat(density))) + scale_x_sqrt(limits = c(0, 1.5e4)) + 
  scale_fill_gradient(name = "", low = "bisque2", high = "darkblue") + 
  facet_wrap(~ sample) + theme_bw() + 
  xlab("Read length") + ylab("Average base quality")

## -------------------------------------------------------------------------- ##
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

png(gsub("\\.rds$", "_nbrreads_sqrt.png", outrds), width = 12, 
    height = 8, unit = "in", res = 400)
print(cowplot::plot_grid(
  ptot,
  cowplot::plot_grid(plengthsqrt, pqual, nrow = 1, 
                     rel_widths = c(1, 1), labels = c("B", "C")),
  ncol = 1, rel_heights = c(1, 1), labels = c("A", ""))
)
dev.off()

png(gsub("\\.rds$", "_nbrreads_linear.png", outrds), width = 12, 
    height = 8, unit = "in", res = 400)
print(cowplot::plot_grid(
  ptot,
  cowplot::plot_grid(plengthlin, pqual, nrow = 1, 
                     rel_widths = c(1, 1), labels = c("B", "C")),
  ncol = 1, rel_heights = c(1, 1), labels = c("A", ""))
)
dev.off()

png(gsub("\\.rds$", "_readlength_vs_quality.png", outrds), width = 11, 
    height = 11, unit = "in", res = 400)
print(plengthqual)
dev.off()

png(gsub("\\.rds$", "_readlength_vs_quality_sqrt.png", outrds), width = 11, 
    height = 11, unit = "in", res = 400)
print(plengthqualsqrt)
dev.off()

png(gsub("\\.rds$", "_readlength_vs_quality_linear.png", outrds), width = 11, 
    height = 11, unit = "in", res = 400)
print(plengthquallin)
dev.off()

## -------------------------------------------------------------------------- ##
saveRDS(NULL, file = outrds)
date()
sessionInfo()

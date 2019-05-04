args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(tximport)
  library(ggplot2)
  library(dplyr)
  library(ggridges)
  library(tibble)
})
source("manuscript_results_figures/Rscripts/remap_sample_names.R")

datasets <- strsplit(datasets, ",")[[1]]
names(datasets) <- datasets
conditions <- strsplit(conditions, ",")[[1]]

print(datasets)
print(conditions)
print(tx2gene)
print(outrds)

tx2gene <- readRDS(tx2gene)

## Read estimated abundances Illumina
files <- list.files("Illumina/salmon31", pattern = "quant.sf", 
                    recursive = TRUE, full.names = TRUE)
names(files) <- basename(dirname(files))
files <- files[names(files) %in% 
                 sample_annotation$sample_orig[sample_annotation$condition %in% conditions]]
salmon <- tximport(files = files, type = "salmon", txOut = TRUE)

salmon$length <- cbind(salmon$length, 
                       txLength = tx2gene$txlength[match(rownames(salmon$length), tx2gene$tx)])
stopifnot(all(rownames(salmon$length) == rownames(salmon$abundance)))

txlengths_illumina <- 
  do.call(dplyr::bind_rows, lapply(colnames(salmon$abundance), function(nm) {
    data.frame(tx_read_length = rep(salmon$length[, "txLength"], 
                                    round(salmon$abundance[, nm] * 10)),
               sample = remap[nm],
               dtype = "Illumina transcript lengths",
               dataset = "Illumina",
               stringsAsFactors = FALSE)
  }))

## Nanopore read lengths
readInfo <- do.call(rbind, lapply(datasets, function(ds) {
  rd <- readRDS(paste0(ds, "/output/", ds, "_nbr_reads.rds"))
  rdf <- rd$fastqs
  rdf <- rdf[!grepl("_orig", names(rdf))]
  rdf <- rdf[names(rdf) %in% 
               sample_annotation$sample_orig[sample_annotation$condition %in% conditions]]
  do.call(dplyr::bind_rows, lapply(names(rdf), function(nm) {
    rdf[[nm]]$reads %>% dplyr::select(read, readLength) %>%
      dplyr::mutate(sample = remap[nm], dataset = remapds[ds]) %>%
      dplyr::left_join(
        rd$genomebams[[nm]]$allAlignments %>%
          dplyr::filter(flag %in% c(0, 16)) %>%
          dplyr::select(read, nbrSupplementaryAlignments) %>%
          dplyr::mutate(alignedGenome = TRUE),
        by = "read"
      ) %>%
      dplyr::select(-read) %>%
      dplyr::mutate(alignedGenome = replace(alignedGenome, is.na(alignedGenome), 
                                            FALSE))
  }))
})) %>%
  dplyr::rename(tx_read_length = readLength) %>%
  dplyr::mutate(dtype = "Nanopore read lengths")

readInfo <- readInfo %>%
  dplyr::mutate(category = "Unaligned") %>%
  dplyr::mutate(category = replace(
    category, alignedGenome & nbrSupplementaryAlignments == 0, 
    "Aligned, without supplementary alignments")) %>%
  dplyr::mutate(category = replace(
    category, alignedGenome & nbrSupplementaryAlignments > 0,
    "Aligned, with supplementary alignment(s)"))

png(gsub("\\.rds$", "_nanoporeread_density_byalignment.png", outrds), 
    width = 7, height = 5, unit = "in", res = 400)
print(ggplot(readInfo %>%
               dplyr::mutate(dataset = factor(dataset, levels = 
                                                ds_order[ds_order %in% dataset])),
             aes(x = tx_read_length, color = category, group = category)) + 
        geom_line(stat = "density", size = 1.25) + theme_bw() + 
        scale_color_manual(values = c("Unaligned" = "lightblue", 
                                      "Aligned, without supplementary alignments" = "#E8601C",
                                      "Aligned, with supplementary alignment(s)" = "#7BAFDE"), 
                           name = "") + 
        scale_x_log10() + xlab("Read length") + 
        theme(legend.position = "bottom") + 
        facet_wrap(~ dataset))
dev.off()

## Plot distribution of transcript lengths from Illumina (weighted by
## abundance), overlay distribution of Nanopore read lengths
plotdf <- dplyr::bind_rows(txlengths_illumina, readInfo %>% dplyr::select(-category)) %>%
  dplyr::mutate(dataset = factor(dataset, levels = ds_order[ds_order %in% dataset]))
## Only aligned reads
plotdfaligned <- dplyr::bind_rows(txlengths_illumina, 
                                  readInfo %>% dplyr::filter(category != "Unaligned") %>%
                                    dplyr::select(-category)) %>%
  dplyr::mutate(dataset = factor(dataset, levels = ds_order[ds_order %in% dataset]))

png(gsub("\\.rds$", "_illuminatx_nanoporeread_density_byds.png", outrds), 
    width = 7, height = 5, unit = "in", res = 400)
dbyds <- ggplot(plotdf,
                aes(x = tx_read_length, color = dataset, group = dataset)) + 
  geom_line(stat = "density", size = 1.25) + theme_bw() + 
  scale_x_log10() + xlab("Transcript/read length") + 
  scale_color_manual(values = ds_colors, name = "")
print(dbyds)
dev.off()

png(gsub("\\.rds$", "_illuminatx_nanoporeread_density_byds_onlyaligned.png", outrds), 
    width = 7, height = 5, unit = "in", res = 400)
dbydsaligned <- ggplot(plotdfaligned,
                       aes(x = tx_read_length, color = dataset, group = dataset)) + 
  geom_line(stat = "density", size = 1.25) + theme_bw() + 
  scale_x_log10() + xlab("Transcript/read length") + 
  scale_color_manual(values = ds_colors, name = "")
print(dbydsaligned)
dev.off()

png(gsub("\\.rds$", "_illuminatx_nanoporeread_density_linear_byds.png", outrds), 
    width = 7, height = 5, unit = "in", res = 400)
dlinearbyds <- ggplot(plotdf,
                      aes(x = tx_read_length, color = dataset, group = dataset)) + 
  geom_line(stat = "density", size = 1.25) + theme_bw() + 
  xlab("Transcript/read length") + xlim(0, 1.5e4) + 
  scale_color_manual(values = ds_colors, name = "")
print(dlinearbyds)
dev.off()

png(gsub("\\.rds$", "_illuminatx_nanoporeread_density_linear_byds_onlyaligned.png", outrds), 
    width = 7, height = 5, unit = "in", res = 400)
dlinearbydsaligned <- ggplot(plotdfaligned,
                             aes(x = tx_read_length, color = dataset, group = dataset)) + 
  geom_line(stat = "density", size = 1.25) + theme_bw() + 
  xlab("Transcript/read length") + xlim(0, 1.5e4) + 
  scale_color_manual(values = ds_colors, name = "")
print(dlinearbydsaligned)
dev.off()

png(gsub("\\.rds$", "_illuminatx_nanoporeread_violin_byds.png", outrds), 
    width = 7, height = 5, unit = "in", res = 400)
dviolinbyds <- ggplot(plotdf,
                      aes(x = dataset, y = tx_read_length,
                          fill = dataset)) + 
  geom_violin(alpha = 0.5) + theme_bw() + 
  scale_y_log10() + ylab("Transcript/read length") + xlab("") + 
  scale_fill_manual(values = ds_colors, name = "") + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(dviolinbyds)
dev.off()

png(gsub("\\.rds$", "_illuminatx_nanoporeread_violin_byds_onlyaligned.png", outrds), 
    width = 7, height = 5, unit = "in", res = 400)
dviolinbydsaligned <- ggplot(plotdfaligned,
                             aes(x = dataset, y = tx_read_length,
                                 fill = dataset)) + 
  geom_violin(alpha = 0.5) + theme_bw() + 
  scale_y_log10() + ylab("Transcript/read length") + xlab("") + 
  scale_fill_manual(values = ds_colors, name = "") + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(dviolinbydsaligned)
dev.off()

png(gsub("\\.rds$", "_illuminatx_nanoporeread_violin_linear_byds.png", outrds), 
    width = 7, height = 5, unit = "in", res = 400)
dviolinlinearbyds <- ggplot(plotdf,
                            aes(x = dataset, y = tx_read_length,
                                fill = dataset)) + 
  geom_violin(alpha = 0.5) + theme_bw() + ylim(0, 1.5e4) + 
  ylab("Transcript/read length") + xlab("") + 
  scale_fill_manual(values = ds_colors, name = "") + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(dviolinlinearbyds)
dev.off()

png(gsub("\\.rds$", "_illuminatx_nanoporeread_violin_linear_byds_onlyaligned.png", outrds), 
    width = 7, height = 5, unit = "in", res = 400)
dviolinlinearbydsaligned <- ggplot(plotdfaligned,
                                   aes(x = dataset, y = tx_read_length,
                                       fill = dataset)) + 
  geom_violin(alpha = 0.5) + theme_bw() + ylim(0, 1.5e4) + 
  ylab("Transcript/read length") + xlab("") + 
  scale_fill_manual(values = ds_colors, name = "") + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(dviolinlinearbydsaligned)
dev.off()

png(gsub("\\.rds$", "_illuminatx_nanoporeread_density_bysample.png", outrds), 
    width = 7, height = 5, unit = "in", res = 400)
ggplot(plotdf,
       aes(x = tx_read_length, color = dataset, group = sample)) + 
  geom_line(stat = "density", size = 1) + theme_bw() + 
  scale_x_log10() + xlab("Transcript/read length") + 
  scale_color_manual(values = ds_colors, name = "")
dev.off()

png(gsub("\\.rds$", "_illuminatx_nanoporeread_ridge_byds.png", outrds), 
    width = 7, height = 7, unit = "in", res = 400)
ggplot(plotdf,
       aes(y = dataset, x = tx_read_length, fill = dataset, 
           color = dataset)) + 
  geom_density_ridges(scale = 1.5) + theme_bw() + 
  scale_x_log10() + xlab("Transcript/read length") + xlab("") + 
  scale_fill_manual(values = ds_colors, name = "") + 
  scale_color_manual(values = ds_colors, name = "")
dev.off()

plots <- list(
  illuminatx_nanoporeread_density_byds = dbyds,
  illuminatx_nanoporeread_density_linear_byds = dlinearbyds,
  illuminatx_nanoporeread_violin_byds = dviolinbyds,
  illuminatx_nanoporeread_violin_linear_byds = dviolinlinearbyds,
  illuminatx_nanoporeread_density_byds_aligned = dbydsaligned,
  illuminatx_nanoporeread_density_linear_byds_aligned = dlinearbydsaligned,
  illuminatx_nanoporeread_violin_byds_aligned = dviolinbydsaligned,
  illuminatx_nanoporeread_violin_linear_byds_aligned = dviolinlinearbydsaligned
)

saveRDS(plots, file = outrds)
date()
sessionInfo()
## Plot summary statistics for data sets

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
    eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tximport)
  library(GenomicAlignments)
  library(cowplot)
})
source("manuscript_results_figures/Rscripts/remap_sample_names.R")

datasets <- strsplit(datasets, ",")[[1]]
names(datasets) <- datasets
conditions <- strsplit(conditions, ",")[[1]]

print(tx2gene)
print(datasets)
print(conditions)
print(txcercc)  ## Plot of transcript coverage, ERCC
print(txcsirv)  ## Plot of transcript coverage, SIRV
print(txcna12878)  ## Plot of transcript coverage, NA12878
print(ilmnnplength)  ## Plot of Illumina/nanopore tx/read length
print(outrds)

tx2gene <- readRDS(tx2gene)

################################################################################
## Number of covered transcripts, covered fraction
################################################################################
dfct_primary <- do.call(dplyr::bind_rows, lapply(datasets, function(ds) {
  rd <- readRDS(paste0(ds, "/output/", ds, "_nbr_reads.rds"))
  rdt <- rd$txomebams_p0.99
  rdt <- rdt[grep("_orig", names(rdt), invert = TRUE)]
  rdt <- rdt[names(rdt) %in% 
               sample_annotation$sample_orig[sample_annotation$condition %in% conditions]]
  do.call(dplyr::bind_rows, lapply(names(rdt), function(nm) {
    rdt[[nm]]$allAlignments %>% dplyr::filter(flag %in% c(0, 16)) %>% 
      dplyr::mutate(sample = remap[nm]) %>% 
      dplyr::mutate(dataset = remapds[ds])
  }))
})) %>%
  dplyr::mutate(gname = tx2gene$gene[match(rname, tx2gene$tx)]) %>%
  dplyr::mutate(dataset = factor(dataset, levels = ds_order[ds_order %in% dataset]))
  
dfct_longest <- do.call(dplyr::bind_rows, lapply(datasets, function(ds) {
  rd <- readRDS(paste0(ds, "/output/", ds, "_nbr_reads.rds"))
  rdt <- rd$txomebams_p0.99
  rdt <- rdt[grep("_orig", names(rdt), invert = TRUE)]
  rdt <- rdt[names(rdt) %in% 
               sample_annotation$sample_orig[sample_annotation$condition %in% conditions]]
  do.call(dplyr::bind_rows, lapply(names(rdt), function(nm) {
    rdt[[nm]]$allAlignments %>% 
      dplyr::select(read, flag, rname, nbrM, nbrS, nbrH, nbrD, nbrI, 
                    txLength, alignedLength, nbrSecondaryAlignments) %>%
      dplyr::filter(flag %in% c(0, 16, 256, 272)) %>% 
      dplyr::group_by(read) %>%
      dplyr::filter(alignedLength/max(alignedLength) > 0.9) %>%
      dplyr::arrange(desc((nbrM + nbrD)/txLength)) %>%
      dplyr::slice(1) %>% 
      dplyr::ungroup() %>% 
      dplyr::mutate(sample = remap[nm]) %>% 
      dplyr::mutate(dataset = remapds[ds])
  }))
})) %>%
  dplyr::mutate(gname = tx2gene$gene[match(rname, tx2gene$tx)]) %>%
  dplyr::mutate(dataset = factor(dataset, levels = ds_order[ds_order %in% dataset]))

## Write to file: for each sample, the fraction of reads aligning to the
## transcriptome that cover at least X% of at least one of the transcripts it
## aligns to
write.table(
  dfct_longest %>% dplyr::mutate(covFraction = (nbrM + nbrD)/txLength) %>% 
    dplyr::group_by(sample, dataset) %>% 
    dplyr::summarize(fracReadsCovAtLeast0.1 = mean(covFraction >= 0.1),
                     fracReadsCovAtLeast0.2 = mean(covFraction >= 0.2),
                     fracReadsCovAtLeast0.3 = mean(covFraction >= 0.3),
                     fracReadsCovAtLeast0.4 = mean(covFraction >= 0.4),
                     fracReadsCovAtLeast0.5 = mean(covFraction >= 0.5),
                     fracReadsCovAtLeast0.6 = mean(covFraction >= 0.6),
                     fracReadsCovAtLeast0.7 = mean(covFraction >= 0.7),
                     fracReadsCovAtLeast0.8 = mean(covFraction >= 0.8),
                     fracReadsCovAtLeast0.9 = mean(covFraction >= 0.9),
                     fracReadsCovAtLeast1.0 = mean(covFraction >= 1.0)),
  file = gsub("\\.rds$", "_fracreads_cov_atleast_xpercent_of_tx_longest.txt", outrds),
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE
)

## Help functions
plotCovFractionOverall <- function(plotdf, maxLength) {
  ggplot(plotdf %>% dplyr::mutate(covFraction = (nbrM + nbrD)/txLength), 
         aes(x = "All transcripts", y = covFraction)) + 
    geom_violin() + scale_y_continuous(limits = c(0, NA)) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    xlab("") + ylab("Fraction of transcript covered per alignment") + 
    stat_summary(data = plotdf %>% 
                   dplyr::summarize(covFraction = length(txLength)), 
                 fun.data = function(x) {return(c(y = 1, label = x))}, 
                 geom = "text", alpha = 1, size = 2.5, vjust = -1)
}

plotCovFractionByLength <- function(plotdf, maxLength) {
  ggplot(plotdf %>% 
           dplyr::mutate(
             txLengthGroup = Hmisc::cut2(txLength, 
                                         cuts = c(0, 500, 1000, 1500, 
                                                  2000, maxLength))) %>%
           dplyr::mutate(covFraction = (nbrM + nbrD)/txLength), 
         aes(x = txLengthGroup, y = covFraction)) + 
    geom_violin() + scale_y_continuous(limits = c(0, NA)) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    xlab("Transcript length") + 
    ylab("Fraction of transcript covered per alignment") + 
    stat_summary(data = plotdf %>% 
                   dplyr::mutate(
                     txLengthGroup = Hmisc::cut2(txLength, 
                                                 cuts = c(0, 500, 1000, 1500, 
                                                          2000, maxLength))) %>%
                   dplyr::group_by(txLengthGroup) %>%
                   dplyr::summarize(covFraction = length(txLength)), 
                 fun.data = function(x) {return(c(y = 1, label = x))}, 
                 geom = "text", alpha = 1, size = 2.5, vjust = -1)
}

plotCovFractionByLengthSample <- function(plotdf, maxLength) {
  ggplot(plotdf, 
         aes(x = Hmisc::cut2(txLength, cuts = c(0, 500, 1000, 1500, 
                                                2000, maxLength)),
             y = (nbrM + nbrD)/txLength)) + 
    facet_wrap(~ sample) + geom_violin() + 
    scale_y_continuous(limits = c(0, NA)) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    xlab("Transcript length") + 
    ylab("Fraction of transcript covered per alignment")
}

plotCovFractionByLengthDataset <- function(plotdf, maxLength) {
  ggplot(plotdf %>% 
           dplyr::mutate(
             txLengthGroup = Hmisc::cut2(txLength, 
                                         cuts = c(0, 500, 1000, 1500, 
                                                  2000, maxLength))) %>%
           dplyr::mutate(covFraction = (nbrM + nbrD)/txLength), 
         aes(x = txLengthGroup, y = covFraction, fill = dataset)) + 
    facet_wrap(~ dataset) + geom_violin(alpha = 0.5) + 
    scale_y_continuous(limits = c(0, NA)) + 
    scale_fill_manual(values = ds_colors, name = "") + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "none") + 
    xlab("Transcript length") + 
    ylab("Fraction of transcript covered per alignment") + 
    stat_summary(data = plotdf %>% 
                   dplyr::mutate(
                     txLengthGroup = Hmisc::cut2(txLength, 
                                                 cuts = c(0, 500, 1000, 1500, 
                                                          2000, maxLength))) %>%
                   dplyr::group_by(dataset, txLengthGroup) %>%
                   dplyr::summarize(covFraction = length(txLength)), 
                 fun.data = function(x) {return(c(y = 0.99, label = x))}, 
                 geom = "text", alpha = 1, size = 2.5, vjust = -1)
}

png(gsub("\\.rds$", "_coverage_fraction_of_transcripts_longest_readswithsecondaryalignments_byds.png", outrds), width = 8, height = 3, unit = "in", res = 400)
print(ggplot(dfct_longest %>% 
               dplyr::mutate(covFraction = (nbrM + nbrD)/txLength) %>%
               dplyr::filter(nbrSecondaryAlignments > 0),
             aes(x = covFraction)) + 
        geom_histogram(bins = 100, aes(fill = dataset)) + theme_bw() + 
        facet_wrap(~ dataset, nrow = 1) + 
        scale_fill_manual(values = ds_colors, name = "") + 
        theme(legend.position = "none") + 
        xlab("Maximal fraction of transcript covered by alignment") + 
        ylab("Count") + ggtitle("Reads with at least one secondary alignment"))
dev.off()

png(gsub("\\.rds$", "_coverage_fraction_of_transcripts_longest_readswithsecondaryalignments_bysample.png", outrds), width = 12, height = 10, unit = "in", res = 400)
print(ggplot(dfct_longest %>% 
               dplyr::mutate(covFraction = (nbrM + nbrD)/txLength) %>%
               dplyr::filter(nbrSecondaryAlignments > 0),
             aes(x = covFraction)) + 
        geom_histogram(bins = 100, aes(fill = dataset)) + theme_bw() + 
        facet_wrap(~ sample) + 
        scale_fill_manual(values = ds_colors, name = "") + 
        theme(legend.position = "none") + 
        xlab("Maximal fraction of transcript covered by alignment") + 
        ylab("Count") + ggtitle("Reads with at least one secondary alignment"))
dev.off()

png(gsub("\\.rds$", "_coverage_fraction_of_transcripts_primary.png", outrds),
    width = 5, height = 5, unit = "in", res = 400)
q1p <- cowplot::plot_grid(
  plotCovFractionOverall(dfct_primary, max(dfct_primary$txLength)),
  plotCovFractionByLength(dfct_primary, max(dfct_primary$txLength)),
  rel_widths = c(2, 5), labels = c("A", "B"))
print(q1p)
dev.off()

png(gsub("\\.rds$", "_coverage_fraction_of_transcripts_longest.png", outrds),
    width = 5, height = 5, unit = "in", res = 400)
q1l <- cowplot::plot_grid(
  plotCovFractionOverall(dfct_longest, max(dfct_longest$txLength)),
  plotCovFractionByLength(dfct_longest, max(dfct_longest$txLength)),
  rel_widths = c(2, 5), labels = c("A", "B"))
print(q1l)
dev.off()

png(gsub("\\.rds$", "_coverage_fraction_of_transcripts_primary_bysample.png", outrds),
    width = 12, height = 10, unit = "in", res = 400)
print(plotCovFractionByLengthSample(dfct_primary, max(dfct_primary$txLength)))
dev.off()

png(gsub("\\.rds$", "_coverage_fraction_of_transcripts_primary_bydataset.png", outrds),
    width = 9, height = 7.5, unit = "in", res = 400)
print(plotCovFractionByLengthDataset(dfct_primary, max(dfct_primary$txLength)))
dev.off()

png(gsub("\\.rds$", "_coverage_fraction_of_transcripts_longest_bysample.png", outrds),
    width = 12, height = 10, unit = "in", res = 400)
print(plotCovFractionByLengthSample(dfct_longest, max(dfct_longest$txLength)))
dev.off()

png(gsub("\\.rds$", "_coverage_fraction_of_transcripts_longest_bydataset.png", outrds),
    width = 9, height = 7.5, unit = "in", res = 400)
pcovlbyds <- plotCovFractionByLengthDataset(dfct_longest, max(dfct_longest$txLength))
print(pcovlbyds)
dev.off()

## Put together with SIRV/ERCC + Illumina comparison
sirvcov <- readRDS(txcsirv)
ercccov <- readRDS(txcercc)
na12878cov <- readRDS(txcna12878)
ilmn <- readRDS(ilmnnplength)

png(gsub("\\.rds$", "_coverage_fraction_summary.png", outrds),
    width = 14, height = 10, unit = "in", res = 400)
print(cowplot::plot_grid(
  plotCovFractionByLengthDataset(dfct_longest, max(dfct_longest$txLength)) + 
    facet_wrap(~ dataset, nrow = 1),
  cowplot::plot_grid(
    na12878cov$covfraclongest + ggtitle("NA12878"),
    sirvcov$covfraclongest + ggtitle("SIRV"),
    ercccov$covfraclongest + ggtitle("ERCC"),
    ilmn$illuminatx_nanoporeread_density_byds + 
      theme(legend.position = c(0.75, 0.77)),
    nrow = 1, rel_widths = c(1.12, 0.8, 0.8, 1.28), labels = c("B", "C", "D", "E"),
    align = "h", axis = "tb"
  ),
  ncol = 1, rel_heights = c(1, 1), labels = c("A", "")
))
dev.off()

png(gsub("\\.rds$", "_coverage_fraction_summary_onlyaligned.png", outrds),
    width = 14, height = 10, unit = "in", res = 400)
print(cowplot::plot_grid(
  plotCovFractionByLengthDataset(dfct_longest, max(dfct_longest$txLength)) + 
    facet_wrap(~ dataset, nrow = 1),
  cowplot::plot_grid(
    na12878cov$covfraclongest + ggtitle("NA12878"),
    sirvcov$covfraclongest + ggtitle("SIRV"),
    ercccov$covfraclongest + ggtitle("ERCC"),
    ilmn$illuminatx_nanoporeread_density_byds_aligned + 
      theme(legend.position = c(0.75, 0.77)),
    nrow = 1, rel_widths = c(1.12, 0.8, 0.8, 1.28), labels = c("B", "C", "D", "E"),
    align = "h", axis = "tb"
  ),
  ncol = 1, rel_heights = c(1, 1), labels = c("A", "")
))
dev.off()

png(gsub("\\.rds$", "_coverage_fraction_summary_linear.png", outrds),
    width = 14, height = 10, unit = "in", res = 400)
print(cowplot::plot_grid(
  plotCovFractionByLengthDataset(dfct_longest, max(dfct_longest$txLength)) + 
    facet_wrap(~ dataset, nrow = 1),
  cowplot::plot_grid(
    na12878cov$covfraclongest + ggtitle("NA12878"),
    sirvcov$covfraclongest + ggtitle("SIRV"),
    ercccov$covfraclongest + ggtitle("ERCC"),
    ilmn$illuminatx_nanoporeread_density_linear_byds + 
      theme(legend.position = c(0.75, 0.77)),
    nrow = 1, rel_widths = c(1.12, 0.8, 0.8, 1.28), labels = c("B", "C", "D", "E"),
    align = "h", axis = "tb"
  ),
  ncol = 1, rel_heights = c(1, 1), labels = c("A", "")
))
dev.off()

png(gsub("\\.rds$", "_coverage_fraction_summary_linear_onlyaligned.png", outrds),
    width = 14, height = 10, unit = "in", res = 400)
print(cowplot::plot_grid(
  plotCovFractionByLengthDataset(dfct_longest, max(dfct_longest$txLength)) + 
    facet_wrap(~ dataset, nrow = 1),
  cowplot::plot_grid(
    na12878cov$covfraclongest + ggtitle("NA12878"),
    sirvcov$covfraclongest + ggtitle("SIRV"),
    ercccov$covfraclongest + ggtitle("ERCC"),
    ilmn$illuminatx_nanoporeread_density_linear_byds_aligned + 
      theme(legend.position = c(0.75, 0.77)),
    nrow = 1, rel_widths = c(1.12, 0.8, 0.8, 1.28), labels = c("B", "C", "D", "E"),
    align = "h", axis = "tb"
  ),
  ncol = 1, rel_heights = c(1, 1), labels = c("A", "")
))
dev.off()

png(gsub("\\.rds$", "_coverage_fraction_summary_violin.png", outrds),
    width = 14, height = 10, unit = "in", res = 400)
print(cowplot::plot_grid(
  plotCovFractionByLengthDataset(dfct_longest, max(dfct_longest$txLength)) + 
    facet_wrap(~ dataset, nrow = 1),
  cowplot::plot_grid(
    na12878cov$covfraclongest + ggtitle("NA12878"),
    sirvcov$covfraclongest + ggtitle("SIRV"),
    ercccov$covfraclongest + ggtitle("ERCC"),
    ilmn$illuminatx_nanoporeread_violin_byds,
    nrow = 1, rel_widths = c(1.12, 0.8, 0.8, 1.28), labels = c("B", "C", "D", "E"),
    align = "h", axis = "tb"
  ),
  ncol = 1, rel_heights = c(1, 1), labels = c("A", "")
))
dev.off()

png(gsub("\\.rds$", "_coverage_fraction_summary_violin_onlyaligned.png", outrds),
    width = 14, height = 10, unit = "in", res = 400)
print(cowplot::plot_grid(
  plotCovFractionByLengthDataset(dfct_longest, max(dfct_longest$txLength)) + 
    facet_wrap(~ dataset, nrow = 1),
  cowplot::plot_grid(
    na12878cov$covfraclongest + ggtitle("NA12878"),
    sirvcov$covfraclongest + ggtitle("SIRV"),
    ercccov$covfraclongest + ggtitle("ERCC"),
    ilmn$illuminatx_nanoporeread_violin_byds_aligned,
    nrow = 1, rel_widths = c(1.12, 0.8, 0.8, 1.28), labels = c("B", "C", "D", "E"),
    align = "h", axis = "tb"
  ),
  ncol = 1, rel_heights = c(1, 1), labels = c("A", "")
))
dev.off()

png(gsub("\\.rds$", "_coverage_fraction_summary_violin_linear.png", outrds),
    width = 14, height = 10, unit = "in", res = 400)
print(cowplot::plot_grid(
  plotCovFractionByLengthDataset(dfct_longest, max(dfct_longest$txLength)) + 
    facet_wrap(~ dataset, nrow = 1),
  cowplot::plot_grid(
    na12878cov$covfraclongest + ggtitle("NA12878"),
    sirvcov$covfraclongest + ggtitle("SIRV"),
    ercccov$covfraclongest + ggtitle("ERCC"),
    ilmn$illuminatx_nanoporeread_violin_linear_byds,
    nrow = 1, rel_widths = c(1.12, 0.8, 0.8, 1.28), labels = c("B", "C", "D", "E"),
    align = "h", axis = "tb"
  ),
  ncol = 1, rel_heights = c(1, 1), labels = c("A", "")
))
dev.off()

png(gsub("\\.rds$", "_coverage_fraction_summary_violin_linear_onlyaligned.png", outrds),
    width = 14, height = 10, unit = "in", res = 400)
print(cowplot::plot_grid(
  plotCovFractionByLengthDataset(dfct_longest, max(dfct_longest$txLength)) + 
    facet_wrap(~ dataset, nrow = 1),
  cowplot::plot_grid(
    na12878cov$covfraclongest + ggtitle("NA12878"),
    sirvcov$covfraclongest + ggtitle("SIRV"),
    ercccov$covfraclongest + ggtitle("ERCC"),
    ilmn$illuminatx_nanoporeread_violin_linear_byds_aligned,
    nrow = 1, rel_widths = c(1.12, 0.8, 0.8, 1.28), labels = c("B", "C", "D", "E"),
    align = "h", axis = "tb"
  ),
  ncol = 1, rel_heights = c(1, 1), labels = c("A", "")
))
dev.off()

palp <- ggplot(dfct_primary, aes(x = txLength, y = nbrM + nbrD)) + 
  facet_wrap(~ sample) + geom_abline(slope = 1, intercept = 0) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  xlab("Transcript length") + ylab("Covered length")

png(gsub("\\.rds$", "_tx_vs_aligned_length_primary_log10_hex.png", outrds),
    width = 10, height = 10, unit = "in", res = 400)
print(palp + geom_hex(bins = 100, aes(fill = stat(density))) +
        scale_fill_gradient(name = "", low = "bisque2", high = "darkblue") + 
        scale_x_log10() + scale_y_log10())
dev.off()

pall <- ggplot(dfct_longest, aes(x = txLength, y = nbrM + nbrD)) + 
  facet_wrap(~ sample) + geom_abline(slope = 1, intercept = 0) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  xlab("Transcript length") + ylab("Covered length")

png(gsub("\\.rds$", "_tx_vs_aligned_length_longest_log10_hex.png", outrds),
    width = 10, height = 10, unit = "in", res = 400)
print(pall + 
        geom_hex(bins = 100, aes(fill = stat(density))) +
        scale_fill_gradient(name = "", low = "bisque2", high = "darkblue") + 
        scale_x_log10() + scale_y_log10())
dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()
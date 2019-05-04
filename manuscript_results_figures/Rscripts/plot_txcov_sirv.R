args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(cowplot)
  library(GenomicAlignments)
})

labels <- strsplit(labels, ",")[[1]]

print(dataset)
print(labels)
print(title)
print(outrds)

################################################################################
## Number of covered transcripts, covered fraction
################################################################################
dfct_primary <- do.call(dplyr::bind_rows, lapply(dataset, function(ds) {
  rd <- readRDS(paste0(ds, "/output/", ds, "_nbr_reads.rds"))
  rdt <- rd$txomebams_p0.99
  do.call(dplyr::bind_rows, lapply(names(rdt), function(nm) {
    rdt[[nm]]$allAlignments %>% dplyr::filter(flag %in% c(0, 16)) %>% 
      dplyr::mutate(sample = nm) %>% 
      dplyr::mutate(dataset = ds)
  }))
}))

dfct_longest <- do.call(dplyr::bind_rows, lapply(dataset, function(ds) {
  rd <- readRDS(paste0(ds, "/output/", ds, "_nbr_reads.rds"))
  rdt <- rd$txomebams_p0.99
  do.call(dplyr::bind_rows, lapply(names(rdt), function(nm) {
    rdt[[nm]]$allAlignments %>% 
      dplyr::select(read, flag, rname, nbrM, nbrS, nbrH, nbrD, nbrI, 
                    txLength, alignedLength) %>%
      dplyr::filter(flag %in% c(0, 16, 256, 272)) %>% 
      dplyr::group_by(read) %>%
      dplyr::filter(alignedLength/max(alignedLength) > 0.9) %>%
      dplyr::arrange(desc((nbrM + nbrD)/txLength)) %>%
      dplyr::slice(1) %>% 
      dplyr::ungroup() %>% 
      dplyr::mutate(sample = nm) %>% 
      dplyr::mutate(dataset = ds)
  }))
}))

plotCovFractionAll <- function(plotdf, maxLength) {
  ggplot(plotdf %>% dplyr::mutate(covFraction = (nbrM + nbrD)/txLength),
         aes(x = "All transcripts", y = covFraction)) + geom_violin() + 
    scale_y_continuous(limits = c(0, NA)) + 
    theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    xlab("") + ylab("Fraction of transcript covered per alignment") + 
    stat_summary(data = plotdf %>% 
                   dplyr::summarize(covFraction = length(txLength)), 
                 fun.data = function(x) {return(c(y = 1, label = x))}, 
                 geom = "text", alpha = 1, size = 2.5, vjust = -1)
}

plotCovFractionStrat <- function(plotdf, maxLength) {
  ggplot(plotdf %>% 
           dplyr::mutate(txLengthGroup = Hmisc::cut2(txLength, 
                                                     cuts = c(0, 500, 1000, 1500, 
                                                              2000, maxLength))) %>%
           dplyr::mutate(covFraction = (nbrM + nbrD)/txLength), 
         aes(x = txLengthGroup, y = covFraction)) + 
    geom_violin() + scale_y_continuous(limits = c(0, NA)) + 
    theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    xlab("Transcript length") + ylab("Fraction of transcript covered per alignment") + 
    stat_summary(data = plotdf %>% 
                   dplyr::mutate(txLengthGroup = Hmisc::cut2(txLength, 
                                                             cuts = c(0, 500, 1000, 1500, 
                                                                      2000, maxLength))) %>%
                   dplyr::group_by(txLengthGroup) %>%
                   dplyr::summarize(covFraction = length(txLength)), 
                 fun.data = function(x) {return(c(y = 1, label = x))}, 
                 geom = "text", alpha = 1, size = 2.5, vjust = -1)
}

png(gsub("\\.rds$", "_coverage_fraction_of_transcripts_primary.png", outrds),
    width = 5, height = 5, unit = "in", res = 400)
plot_title <- cowplot::ggdraw() + cowplot::draw_label(title, fontface = 'bold')
p0 <- cowplot::plot_grid(
  plotCovFractionAll(dfct_primary, max(dfct_primary$txLength)),
  plotCovFractionStrat(dfct_primary, max(dfct_primary$txLength)),
  rel_widths = c(2, 5), labels = labels
)
print(cowplot::plot_grid(plot_title, p0, rel_heights = c(0.075, 1), labels = "", ncol = 1))
dev.off()

png(gsub("\\.rds$", "_coverage_fraction_of_transcripts_longest.png", outrds),
    width = 5, height = 5, unit = "in", res = 400)
plot_title <- cowplot::ggdraw() + cowplot::draw_label(title, fontface = 'bold')
p0 <- cowplot::plot_grid(
  plotCovFractionAll(dfct_longest, max(dfct_longest$txLength)),
  plotCovFractionStrat(dfct_longest, max(dfct_longest$txLength)),
  rel_widths = c(2, 5), labels = labels
)
print(cowplot::plot_grid(plot_title, p0, rel_heights = c(0.075, 1), labels = "", ncol = 1))
dev.off()

saveRDS(list(covfraclongest = plotCovFractionStrat(dfct_longest, 
                                                   max(dfct_longest$txLength))), 
        file = outrds)
date()
sessionInfo()

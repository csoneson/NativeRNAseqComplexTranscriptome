args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tximport)
})
source("manuscript_results_figures/Rscripts/remap_sample_names.R")

datasets <- strsplit(datasets, ",")[[1]]
names(datasets) <- datasets
conditions <- strsplit(conditions, ",")[[1]]

print(datasets)
print(conditions)
print(outrds)

## Read data
dfprimsupp <- do.call(dplyr::bind_rows, lapply(datasets, function(ds) {
  readRDS(paste0(ds, "/output/", ds, 
                 "_primary_supplementary_alignments_distances.rds")) %>% 
    dplyr::mutate(sample = remap[sample]) %>% 
    dplyr::mutate(dataset = remapds[ds])
})) %>%
  dplyr::left_join(sample_annotation %>% dplyr::select(sample_remap, condition),
                   by = c("sample" = "sample_remap")) %>%
  dplyr::filter(condition %in% conditions) %>%
  dplyr::select(-condition) %>%
  dplyr::mutate(dataset = factor(dataset, levels = ds_order[ds_order %in% dataset]))

g0 <- ggplot(dfprimsupp %>% 
               dplyr::mutate(sample = removeDatasetFromSample(sample, dataset)) %>%
               dplyr::mutate(dist0 = "") %>%
               dplyr::mutate(dist0 = replace(dist0, distn == 0, ", overlapping"),
                             dist0 = replace(dist0, distn > 0, ", non-overlapping")) %>% 
               dplyr::group_by(sample, dataset, strands, dist0) %>% dplyr::tally() %>%
               dplyr::mutate(strands_dist0 = paste0(strands, dist0)),
             aes(x = sample, y = n, fill = strands_dist0)) + 
  facet_grid(~ dataset, scales = "free_x", space = "free_x") + 
  theme_bw() + 
  theme(legend.position = "bottom",
        strip.text = element_text(size = 8)) + 
  xlab("") + 
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) + 
  scale_fill_manual(values = c(
    "different chromosomes" = "#E8601C",
    "same chromosome, different strand, non-overlapping" = "#7BAFDE", 
    "same chromosome, different strand, overlapping" = "#90C987", 
    "same chromosome, same strand, non-overlapping" = "#777777", 
    "same chromosome, same strand, overlapping" = "#B17BA6"
  ), name = "")

## By dataset
g0ds <- ggplot(dfprimsupp %>% 
               dplyr::mutate(dist0 = "") %>%
               dplyr::mutate(dist0 = replace(dist0, distn == 0, ", overlapping"),
                             dist0 = replace(dist0, distn > 0, ", non-overlapping")) %>% 
               dplyr::group_by(dataset, strands, dist0) %>% dplyr::tally() %>%
               dplyr::mutate(strands_dist0 = paste0(strands, dist0)),
             aes(x = dataset, y = n, fill = strands_dist0)) + 
  theme_bw() + 
  theme(legend.position = "bottom",
        strip.text = element_text(size = 8),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  xlab("") + 
  guides(fill = guide_legend(nrow = 5, byrow = TRUE)) + 
  scale_fill_manual(values = c(
    "different chromosomes" = "#E8601C",
    "same chromosome, different strand, non-overlapping" = "#7BAFDE", 
    "same chromosome, different strand, overlapping" = "#90C987", 
    "same chromosome, same strand, non-overlapping" = "#777777", 
    "same chromosome, same strand, overlapping" = "#B17BA6"
  ), name = "")


png(gsub("\\.rds$", "_nbr_primary_supplementary_pairs.png", outrds),
    width = 10, height = 6, unit = "in", res = 400)
print(g0 + geom_bar(stat = "identity") + 
        scale_y_continuous(expand = c(0, 0, 0.05, 0)) + 
        ylab("Number of primary/supplementary\nalignment pairs"))
dev.off()

png(gsub("\\.rds$", "_frac_primary_supplementary_pairs.png", outrds),
    width = 10, height = 6, unit = "in", res = 400)
print(g0 + geom_bar(stat = "identity", position = "fill") + 
        scale_y_continuous(expand = c(0, 0, 0.05, 0)) + 
        ylab("Fraction of primary/supplementary\nalignment pairs"))
dev.off()

plots <- list(
  primary_supplementary_fraction = 
    g0 + geom_bar(stat = "identity", position = "fill") + 
    scale_y_continuous(expand = c(0, 0, 0.05, 0)) + 
    ylab("Fraction of primary/supplementary\nalignment pairs"),
  primary_supplementary_fraction_byds = 
    g0ds + geom_bar(stat = "identity", position = "fill") + 
    scale_y_continuous(expand = c(0, 0, 0.05, 0)) + 
    ylab("Fraction of primary/supplementary\nalignment pairs")
)

saveRDS(plots, file = outrds)
date()
sessionInfo()

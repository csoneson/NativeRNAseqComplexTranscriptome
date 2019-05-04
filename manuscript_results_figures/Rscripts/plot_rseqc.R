args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(tidyr)
})
source("manuscript_results_figures/Rscripts/remap_sample_names.R")
source("RSeQC/RSeQC_geneBody_coverage.geneBodyCoverage.r")

datasets <- strsplit(datasets, ",")[[1]]
datasets <- c(datasets, "Illumina")
conditions <- strsplit(conditions, ",")[[1]]

print(datasets)
print(conditions)
print(outrds)

df <- data.frame(sample = gsub("0918\\.A_", "0918\\.A-", gsub("^V", "", rowLabel)), 
                 data.frame(data_matrix) %>% setNames(paste0("x", seq_len(100))),
                 stringsAsFactors = FALSE) %>%
  tidyr::gather(key = "x", value = "y", -sample) %>%
  dplyr::mutate(x = as.numeric(gsub("x", "", x)),
                sample = remap[gsub("_minimap_genome_s|_Aligned.sortedByCoord.out",
                                    "", sample)]) %>%
  dplyr::left_join(sample_annotation %>% 
                     dplyr::select(sample_remap, condition, dataset), 
                   by = c("sample" = "sample_remap")) %>%
  dplyr::filter(dataset %in% remapds[datasets] & condition %in% conditions) %>% 
  dplyr::mutate(
    dataset = factor(dataset, levels = ds_order[ds_order %in% dataset])
  )

png(gsub("\\.rds$", "_gene_body_coverage.png", outrds), width = 6, height = 3.5,
    unit = "in", res = 400)
ggplot(df, aes(x = x, y = y, group = sample, color = dataset)) + 
  geom_line() + theme_bw() + 
  scale_color_manual(values = ds_colors, name = "") + 
  guides(color = guide_legend(override.aes = list(size = 2))) + 
  xlab("Gene body percentile (5' -> 3')") + ylab("Coverage")
dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()

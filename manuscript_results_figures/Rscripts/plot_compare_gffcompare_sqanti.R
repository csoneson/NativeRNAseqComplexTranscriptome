args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(rtracklayer)
  library(dplyr)
  library(GenomicFeatures)
  library(ggplot2)
  library(grDevices)
  library(cowplot)
  library(grid)
  library(pheatmap)
  library(viridis)
})
source("manuscript_results_figures/Rscripts/remap_sample_names.R")

datasets <- strsplit(datasets, ",")[[1]]
names(datasets) <- datasets
conditions <- strsplit(conditions, ",")[[1]]

if (all(c("wt", "srpk") %in% conditions)) {
  flaircond <- "_all"
} else if ("wt" %in% conditions) {
  flaircond <- "_WT"
} else {
  stop("Unknown conditions")
}

print(datasets)
print(conditions)
print(flaircond)
print(outrds)

## Read class codes from gffcompare
dfcc <- do.call(dplyr::bind_rows, lapply(datasets, function(ds) {
  tmap <- read.delim(paste0(
    ds, "/flair_round2/", ds, flaircond, "/", ds, flaircond, 
    "_minimap_genome_s_primary_flair_collapse.isoforms.gffcompare.", 
    ds, flaircond, "_minimap_genome_s_primary_flair_collapse.isoforms.gtf.tmap"), 
    header = TRUE, as.is = TRUE) %>%
    dplyr::select(qry_id, class_code, ref_id) %>%
    dplyr::mutate(dataset = remapds[ds]) %>%
    dplyr::mutate(class_code = as.character(class_code),
                  ref_id = as.character(ref_id))
  if (ds %in% c("RNA001", "DCS108", "pilot")) {
    tmap2 <- read.delim(paste0(
      ds, "/flair_round2_ilmnjunc/", ds, flaircond, "/", ds, flaircond, 
      "_minimap_genome_s_primary_flair_collapse_ilmnjunc.isoforms.gffcompare.", 
      ds, flaircond, "_minimap_genome_s_primary_flair_collapse_ilmnjunc.isoforms.",
      "gtf.tmap"), 
      header = TRUE, as.is = TRUE) %>%
      dplyr::select(qry_id, class_code, ref_id) %>%
      dplyr::mutate(dataset = paste0(remapds[ds], "_ILMNjunc")) %>%
      dplyr::mutate(class_code = as.character(class_code),
                    ref_id = as.character(ref_id))
    tmap <- dplyr::bind_rows(tmap, tmap2)
  }
  tmap
})) %>%
  dplyr::mutate(
    class_code = replace(class_code, class_code == "=", 
                         "complete, exact match of intron chain (=)"),
    class_code = replace(class_code, class_code == "c", 
                         "contained in reference (intron compatible) (c)"),
    class_code = replace(class_code, class_code == "k", 
                         "containment of reference (k)"),
    class_code = replace(class_code, class_code == "m",
                         "retained intron(s), full intron chain overlap (m)"),
    class_code = replace(class_code, class_code == "n", 
                         "retained intron(s), partial or no intron chain match (n)"),
    class_code = replace(class_code, class_code == "j", 
                         "multi-exon with at least one junction match (j)"),
    class_code = replace(class_code, class_code == "e", 
                         "single exon transfrag partially covering an intron (e)"),
    class_code = replace(class_code, class_code == "o", 
                         "other same strand overlap with reference exons (o)"),
    class_code = replace(class_code, class_code == "s", 
                         "intron matching on opposite strand (s)"),
    class_code = replace(class_code, class_code == "x", 
                         "exonic overlap on opposite strand (x)"),
    class_code = replace(class_code, class_code == "i", 
                         "fully contained within a reference intron (i)"),
    class_code = replace(class_code, class_code == "y", 
                         "contains a reference within its introns (y)"),
    class_code = replace(class_code, class_code == "p", 
                         "possible polymerase run-on (no overlap) (p)"),
    class_code = replace(class_code, class_code == "r", "repeat (r)"),
    class_code = replace(class_code, class_code == "u", "other (u)")) %>%
  dplyr::mutate(class_code = factor(class_code, levels = c(
    "other (u)",
    "repeat (r)",
    "possible polymerase run-on (no overlap) (p)",
    "contains a reference within its introns (y)",
    "fully contained within a reference intron (i)",
    "exonic overlap on opposite strand (x)",
    "intron matching on opposite strand (s)",
    "other same strand overlap with reference exons (o)",
    "single exon transfrag partially covering an intron (e)",
    "multi-exon with at least one junction match (j)",
    "retained intron(s), partial or no intron chain match (n)",
    "retained intron(s), full intron chain overlap (m)",
    "containment of reference (k)",
    "contained in reference (intron compatible) (c)",
    "complete, exact match of intron chain (=)"
  )))

dfsc <- do.call(dplyr::bind_rows, lapply(datasets, function(ds) {
  tmap <- read.delim(paste0(
    ds, "/sqanti_round2/", ds, flaircond, "/", ds, flaircond, 
    "_classification.txt"), 
    header = TRUE, as.is = TRUE) %>%
    dplyr::select(isoform, structural_category, associated_transcript) %>%
    dplyr::mutate(dataset = remapds[ds]) %>%
    dplyr::mutate(structural_category = as.character(structural_category),
                  associated_transcript = as.character(associated_transcript))
  if (ds %in% c("RNA001", "DCS108", "pilot")) {
    tmap2 <- read.delim(paste0(
      ds, "/sqanti_round2_ilmnjunc/", ds, flaircond, "/", ds, flaircond, 
      "_ilmnjunc_classification.txt"), 
      header = TRUE, as.is = TRUE) %>%
      dplyr::select(isoform, structural_category, associated_transcript) %>%
      dplyr::mutate(dataset = paste0(remapds[ds], "_ILMNjunc")) %>%
      dplyr::mutate(structural_category = as.character(structural_category),
                    associated_transcript = as.character(associated_transcript))
    tmap <- dplyr::bind_rows(tmap, tmap2)
  }
  tmap
})) %>%
  dplyr::mutate(structural_category = factor(structural_category, levels = c(
    "intergenic",
    "genic_intron",
    "antisense",
    "novel_not_in_catalog",
    "novel_in_catalog",
    "genic",
    "fusion",
    "incomplete-splice_match",
    "full-splice_match"
  )))

dfboth <- dplyr::full_join(dfcc, dfsc, by = c("qry_id" = "isoform", "dataset"))
dfbothsub <- dfboth %>% dplyr::filter(dataset == "ONT-RNA001-HAP_ILMNjunc")

tb <- table(dfboth$class_code, dfboth$structural_category)
tbsub <- table(dfbothsub$class_code, dfbothsub$structural_category)

sqrt_breaks <- function(xs, n = 10) {
  breaks <- (seq(0, (max(xs) + 1)^(1/3), length.out = n))^3
  breaks[!duplicated(breaks)]
}
mat_breaks <- sqrt_breaks(tb, n = 101)
mat_breaks_sub <- sqrt_breaks(tbsub, n = 101)

png(gsub("rds$", "png", outrds), height = 8, width = 8, unit = "in", res = 400)
pheatmap::pheatmap(tb, cluster_rows = FALSE, cluster_cols = FALSE,
                   color = inferno(length(mat_breaks) - 1),
                   breaks = mat_breaks)
dev.off()

png(gsub("\\.rds$", "_sub.png", outrds), height = 8, width = 8, unit = "in", res = 400)
pheatmap::pheatmap(tbsub, cluster_rows = FALSE, cluster_cols = FALSE,
                   color = inferno(length(mat_breaks_sub) - 1),
                   breaks = mat_breaks_sub, main = "ONT-RNA001-HAP_ILMNjunc")
dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()

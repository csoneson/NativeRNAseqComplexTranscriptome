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
})
source("manuscript_results_figures/Rscripts/remap_sample_names.R")

datasets <- strsplit(datasets, ",")[[1]]
datasets <- c(datasets, "Illumina")
names(datasets) <- datasets
conditions <- strsplit(conditions, ",")[[1]]

if (all(c("wt", "srpk") %in% conditions)) {
  sqanticond <- "_all"
} else if ("wt" %in% conditions) {
  sqanticond <- "_WT"
} else {
  stop("Unknown conditions")
}

print(datasets)
print(conditions)
print(sqanticond)
print(outrds)

muted <- c("#DC050C", "#E8601C", "#7BAFDE", "#1965B0", "#B17BA6",
           "#882E72", "#F1932D", "#F6C141", "#F7EE55", "#4EB265",
           "#CAEDAB", "#777777")
colfun <- grDevices::colorRampPalette(muted)

## Read Illumina abundance estimates
ilmn_tx_abundances <- 
  readRDS("Illumina/output/Illumina_all_abundances.rds")$tx_abundances
ilmn_samples <- sample_annotation %>% 
  dplyr::filter(dataset == "Illumina" & condition %in% conditions) %>% 
  dplyr::pull("sample_orig")
print(ilmn_samples)
ave_ilmn_tx_tpm <- 
  data.frame(transcript_id = rownames(ilmn_tx_abundances), 
             ave_tx_tpm = rowMeans(
               ilmn_tx_abundances[, paste0(ilmn_samples, "__tpm__salmon")]
             ),
             stringsAsFactors = FALSE)

## ========================================================================== ##
## Plot overlap of detected transcripts
df0 <- do.call(dplyr::bind_rows, lapply(setdiff(datasets, "Illumina"), function(ds) {
  dplyr::bind_rows(
    read.delim(paste0(
      ds, "/sqanti_round2_ilmnjunc/", ds, sqanticond, "/", ds, sqanticond, 
      "_ilmnjunc_classification.txt"), 
      header = TRUE, as.is = TRUE) %>%
      dplyr::filter(structural_category %in% c("full-splice_match", 
                                               "incomplete-splice_match")) %>%
      dplyr::select(associated_transcript) %>%
      dplyr::rename(transcript_id = associated_transcript) %>%
      dplyr::distinct() %>%
      dplyr::mutate(dataset = paste0(remapds[ds], "_ILMNjunc"),
                    identified = TRUE),
    read.delim(paste0(
      ds, "/sqanti_round2/", ds, sqanticond, "/", ds, sqanticond, 
      "_classification.txt"), 
      header = TRUE, as.is = TRUE) %>%
      dplyr::filter(structural_category %in% c("full-splice_match", 
                                               "incomplete-splice_match")) %>%
      dplyr::select(associated_transcript) %>%
      dplyr::rename(transcript_id = associated_transcript) %>%
      dplyr::distinct() %>%
      dplyr::mutate(dataset = remapds[ds],
                    identified = TRUE)
  )
})) %>%
  dplyr::filter(dataset != "ONT-RNA001-HEK_ILMNjunc")  ## don't use HAP junctions for HEK
dettx <- unique(df0$transcript_id)
detmat <- matrix(0, nrow = length(dettx), ncol = length(unique(df0$dataset)))
rownames(detmat) <- dettx
colnames(detmat) <- unique(df0$dataset)
detmat[as.matrix(df0 %>% dplyr::select(transcript_id, dataset) %>%
                   dplyr::mutate(transcript_id = as.character(transcript_id),
                                 dataset = as.character(dataset)))] <- 1
ordr <- ds_order[ds_order %in% colnames(detmat)]
ordr <- rep(ordr, each = 2)
suffx <- rep(c("", "_ILMNjunc"), length(ordr)/2)
ordr <- paste0(ordr, suffx)
ordr <- setdiff(ordr, "ONT-RNA001-HEK_ILMNjunc")
stopifnot(ordr %in% colnames(detmat),
          colnames(detmat) %in% ordr)
detmat <- detmat[, ordr]
png(gsub("\\.rds$", "_shared_transcripts.png", outrds), 
    width = 10, height = 6, unit = "in", res = 400)
UpSetR::upset(data.frame(detmat, check.names = FALSE), order.by = "freq",
              decreasing = TRUE, keep.order = TRUE, nsets = length(ordr), 
              mainbar.y.label = "Number of shared annotated transcripts",
              sets.x.label = "Number of identified\nannotated transcripts",
              sets.bar.color = ds_colors[gsub("_ILMNjunc", "", colnames(detmat))],
              sets = colnames(detmat), mb.ratio = c(0.55, 0.45))
grid::grid.edit("arrange", name = "arrange2")
upsetplot <- grid::grid.grab()
dev.off()

## ========================================================================== ##
## Plot abundance of identified/unidentified transcripts 
df1 <- df0 %>% 
  dplyr::full_join(
    do.call(dplyr::bind_rows, lapply(setdiff(datasets, "Illumina"), function(ds) {
      dplyr::bind_rows(ave_ilmn_tx_tpm %>% dplyr::mutate(dataset = remapds[ds]),
                       ave_ilmn_tx_tpm %>% dplyr::mutate(dataset = paste0(remapds[ds],
                                                                          "_ILMNjunc"))) %>%
        dplyr::filter(dataset != "ONT-RNA001-HEK_ILMNjunc")
    })), 
    by = c("transcript_id", "dataset")
  ) %>%
  dplyr::mutate(identified = replace(identified, is.na(identified), FALSE))

png(gsub("\\.rds$", "_illumina_abundance_by_identification.png", outrds),
    width = 8, height = 4, unit = "in", res = 400)
ggplot(df1, aes(x = identified, y = ave_tx_tpm + 1)) + 
  geom_boxplot() + facet_wrap(~ dataset, nrow = 1) + theme_bw() + 
  scale_y_log10() + xlab("Transcript identified by FLAIR") + 
  ylab("Average abundance across Illumina samples (TPM + 1)")
dev.off()

## Merge with upset plot
df2 <- df1 %>% dplyr::filter(grepl("_ILMNjunc$|RNA001-HEK", dataset)) %>%
  dplyr::select(-dataset) %>%
  dplyr::group_by(transcript_id, ave_tx_tpm) %>%
  dplyr::summarize(identified = any(identified))
p1 <- ggplot(df2, aes(x = identified, y = ave_tx_tpm + 1)) + 
  geom_boxplot() + theme_bw() + 
  scale_y_log10() + xlab("Transcript identified by FLAIR") + 
  ylab("Average abundance across Illumina samples (TPM + 1)")
png(gsub("\\.rds$", "_illumina_abundance_plus_shared_tx.png", outrds),
    width = 14, height = 6, unit = "in", res = 400)
print(cowplot::plot_grid(
  upsetplot,
  p1, 
  rel_widths = c(5, 1.5), nrow = 1, labels = c("A", "B")
))
dev.off()

## ========================================================================== ##
## Get number of exons and read length per flair transcript and plot for
## different structural categories
ds_order_both <- rep(ds_order, each = 2)
suffx <- rep(c("", "_ILMNjunc"), length(ds_order_both)/2)
ds_order_both <- paste0(ds_order_both, suffx)
nbrexons <- do.call(dplyr::bind_rows, lapply(datasets, function(ds) {
  tmap <- read.delim(paste0(
    ds, "/sqanti", ifelse(ds == "Illumina", "/", "_round2/"),
    ds, sqanticond, "/", ds, sqanticond, 
    "_classification.txt"), 
    header = TRUE, as.is = TRUE) %>%
    dplyr::mutate(dataset = remapds[ds])
  if (ds %in% c("RNA001", "DCS108", "pilot")) {
    tmap2 <- read.delim(paste0(
      ds, "/sqanti_round2_ilmnjunc/", ds, sqanticond, "/", ds, sqanticond, 
      "_ilmnjunc_classification.txt"), 
      header = TRUE, as.is = TRUE) %>%
      dplyr::mutate(dataset = paste0(remapds[ds], "_ILMNjunc"))
    tmap <- dplyr::bind_rows(tmap, tmap2)
  }
  tmap
})) %>%
  dplyr::mutate(exons = replace(exons, exons > 2, ">2")) %>%
  dplyr::mutate(dataset = factor(dataset, levels = ds_order_both[ds_order_both %in% dataset]))

## ========================================================================== ##
## Plot distribution of structural categories
dfcc <- do.call(dplyr::bind_rows, lapply(datasets, function(ds) {
  tmap <- read.delim(paste0(
    ds, "/sqanti", ifelse(ds == "Illumina", "/", "_round2/"),
    ds, sqanticond, "/", ds, sqanticond, 
    "_classification.txt"), 
    header = TRUE, as.is = TRUE)
  t1 <- as.data.frame(table(tmap$structural_category)) %>% 
    dplyr::mutate(dataset = remapds[ds]) %>%
    dplyr::rename(structural_category = Var1) %>%
    dplyr::mutate(structural_category = as.character(structural_category))
  if (ds %in% c("RNA001", "DCS108", "pilot")) {
    tmap2 <- read.delim(paste0(
      ds, "/sqanti_round2_ilmnjunc/", ds, sqanticond, "/", ds, sqanticond, 
      "_ilmnjunc_classification.txt"), 
      header = TRUE, as.is = TRUE)
    t1 <- dplyr::bind_rows(
      t1, as.data.frame(table(tmap2$structural_category)) %>% 
      dplyr::mutate(dataset = paste0(remapds[ds], "_ILMNjunc")) %>%
      dplyr::rename(structural_category = Var1) %>%
      dplyr::mutate(structural_category = as.character(structural_category)))
  }
  t1
})) %>%
  dplyr::mutate(structural_category = factor(structural_category, levels = c(
    "novel_not_in_catalog",
    "novel_in_catalog",
    "antisense",
    "intergenic",
    "genic",
    "genic_intron",
    "fusion",
    "incomplete-splice_match",
    "full-splice_match"
  ))) %>%
  dplyr::mutate(dataset = factor(dataset, levels = ds_order_both[ds_order_both %in% dataset]))

p0 <- ggplot(dfcc, aes(x = dataset, y = Freq)) + 
  geom_bar(stat = "identity", position = "fill", 
           aes(fill = structural_category)) + theme_bw() + 
  scale_fill_manual(values = structure(colfun(15), 
                                       names = levels(dfcc$structural_category)),
                    name = "Overlap type (structural category)") + 
  ylab("Fraction of transcripts") + xlab("") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  scale_y_continuous(expand = c(0, 0, 0.05, 0)) + 
  geom_text(data = dfcc %>% dplyr::group_by(dataset) %>%
              dplyr::summarize(nbrTx = sum(Freq)),
            aes(x = dataset, y = 1, label = nbrTx), 
            vjust = 0, nudge_y = 0.02, size = 3)

png(gsub("\\.rds$", "_nbr_exons_read_length_by_structuralcategory_full.png", outrds), 
    width = 11, height = 16, unit = "in", res = 400)
p1 <- ggplot(nbrexons %>%
               dplyr::group_by(dataset, structural_category, exons) %>%
               dplyr::tally() %>%
               dplyr::ungroup() %>%
               dplyr::mutate(exons = factor(exons, 
                                            levels = c(">2", "2", "1"))), 
             aes(x = structural_category, y = n, fill = exons)) + 
  geom_bar(stat = "identity", position = "fill") + theme_bw() + 
  facet_wrap(~ dataset, nrow = 2) +
  xlab("Structural category") + ylab("Fraction of transcripts") + 
  scale_fill_manual(values = c("1" = "red", "2" = "blue", ">2" = "grey"),
                    name = "Number of\nexons") + 
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  scale_y_continuous(expand = c(0, 0, 0.05, 0))
p2 <- ggplot(nbrexons, 
             aes(x = structural_category, y = length)) + 
  geom_boxplot() + theme_bw() + 
  facet_wrap(~ dataset, nrow = 2) +
  scale_y_log10() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  xlab("Structural category") + ylab("Inferred transcript length")
cowplot::plot_grid(p0, p1, p2, ncol = 1, labels = c("A", "B", "C"), 
                   rel_heights = c(1, 1.1, 1), align = "v", axis = "l")
dev.off()

png(gsub("\\.rds$", "_nbr_exons_read_length_by_structuralcategory.png", outrds), 
    width = 7, height = 12, unit = "in", res = 400)
p1 <- ggplot(nbrexons %>%
               dplyr::filter(dataset %in% c("ONT-RNA001-HAP_ILMNjunc", "Illumina")) %>%
               dplyr::group_by(dataset, structural_category, exons) %>%
               dplyr::tally() %>%
               dplyr::ungroup() %>%
               dplyr::mutate(exons = factor(exons, 
                                            levels = c(">2", "2", "1"))), 
             aes(x = structural_category, y = n, fill = exons)) + 
  geom_bar(stat = "identity", position = "fill") + theme_bw() + 
  facet_wrap(~ dataset, nrow = 1) +
  xlab("Structural category") + ylab("Fraction of transcripts") + 
  scale_fill_manual(values = c("1" = "red", "2" = "blue", ">2" = "grey"),
                    name = "Number of\nexons") + 
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  scale_y_continuous(expand = c(0, 0, 0.05, 0))
p2 <- ggplot(nbrexons %>%
               dplyr::filter(dataset %in% c("ONT-RNA001-HAP_ILMNjunc", "Illumina")), 
             aes(x = structural_category, y = length)) + 
  geom_boxplot() + theme_bw() + 
  facet_wrap(~ dataset, nrow = 1) +
  scale_y_log10() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  xlab("Structural category") + ylab("Inferred transcript length")
cowplot::plot_grid(p0, p1, p2,
                   ncol = 1, labels = c("A", "B", "C"), 
                   rel_heights = c(1.35, 1.28, 1.1), align = "v", axis = "l")
dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()
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
ave_ilmn_tx_tpm <- 
  data.frame(transcript_id = rownames(ilmn_tx_abundances), 
             ave_tx_tpm = rowMeans(
               ilmn_tx_abundances[, paste0(ilmn_samples, "__tpm__salmon")]
             ),
             stringsAsFactors = FALSE)

## ========================================================================== ##
## Plot overlap of detected transcripts
df0 <- do.call(dplyr::bind_rows, lapply(setdiff(datasets, "Illumina"), function(ds) {
  read.delim(paste0(
    ds, "/sqanti/", ds, sqanticond, "/", ds, sqanticond, 
    "_classification.txt"), 
    header = TRUE, as.is = TRUE) %>%
    dplyr::filter(structural_category %in% c("full-splice_match", 
                                             "incomplete-splice_match")) %>%
    dplyr::select(associated_transcript) %>%
    dplyr::rename(transcript_id = associated_transcript) %>%
    dplyr::distinct() %>%
    dplyr::mutate(dataset = remapds[ds],
                  identified = TRUE)
}))
dettx <- unique(df0$transcript_id)
detmat <- matrix(0, nrow = length(dettx), ncol = length(unique(df0$dataset)))
rownames(detmat) <- dettx
colnames(detmat) <- unique(df0$dataset)
detmat[as.matrix(df0 %>% dplyr::select(transcript_id, dataset))] <- 1
png(gsub("\\.rds$", "_shared_transcripts.png", outrds), 
    width = 10, height = 6, unit = "in", res = 400)
UpSetR::upset(data.frame(detmat, check.names = FALSE), order.by = "freq")
grid::grid.edit("arrange", name = "arrange2")
upsetplot <- grid::grid.grab()
dev.off()

## ========================================================================== ##
## Plot abundance of identified/unidentified transcripts 
df1 <- df0 %>% 
  dplyr::full_join(
    do.call(dplyr::bind_rows, lapply(setdiff(datasets, "Illumina"), function(ds) {
      ave_ilmn_tx_tpm %>% dplyr::mutate(dataset = remapds[ds])
    })), 
    by = c("transcript_id", "dataset")
  ) %>%
  dplyr::mutate(identified = replace(identified, is.na(identified), FALSE))

png(gsub("\\.rds$", "_illumina_abundance_by_identification.png", outrds),
    width = 8, height = 4, unit = "in", res = 400)
ggplot(df1, aes(x = identified, y = ave_tx_tpm + 1)) + 
  geom_boxplot() + facet_wrap(~ dataset, nrow = 1) + theme_bw() + 
  scale_y_log10() + xlab("Transcript identified by flair") + 
  ylab("Average abundance across Illumina samples (TPM + 1)")
dev.off()

## Merge with upset plot
df2 <- df1 %>% dplyr::select(-dataset) %>%
  dplyr::group_by(transcript_id, ave_tx_tpm) %>%
  dplyr::summarize(identified = any(identified))
p1 <- ggplot(df2, aes(x = identified, y = ave_tx_tpm + 1)) + 
  geom_boxplot() + theme_bw() + 
  scale_y_log10() + xlab("Transcript identified by flair") + 
  ylab("Average abundance across Illumina samples (TPM + 1)")
png(gsub("\\.rds$", "_illumina_abundance_plus_shared_tx.png", outrds),
    width = 12, height = 6, unit = "in", res = 400)
print(cowplot::plot_grid(
  upsetplot,
  p1, 
  rel_widths = c(4, 1.5), nrow = 1, labels = c("A", "B")
))
dev.off()

## ========================================================================== ##
## Get number of exons and read length per flair transcript and plot for
## different structural categories
nbrexons <- do.call(dplyr::bind_rows, lapply(datasets, function(ds) {
  read.delim(paste0(
    ds, "/sqanti/", ds, sqanticond, "/", ds, sqanticond, "_classification.txt"), 
    header = TRUE, as.is = TRUE) %>%
    dplyr::mutate(dataset = remapds[ds])
})) %>%
  dplyr::mutate(exons = replace(exons, exons > 2, ">2"))

## ========================================================================== ##
## Plot distribution of structural categories
dfcc <- do.call(dplyr::bind_rows, lapply(datasets, function(ds) {
  tmap <- read.delim(paste0(
    ds, "/sqanti/", ds, sqanticond, "/", ds, sqanticond, "_classification.txt"), 
    header = TRUE, as.is = TRUE)
  as.data.frame(table(tmap$structural_category)) %>% 
    dplyr::mutate(dataset = remapds[ds]) %>%
    dplyr::rename(structural_category = Var1) %>%
    dplyr::mutate(structural_category = as.character(structural_category))
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
  )))

png(gsub("\\.rds$", "_nbr_exons_read_length_by_structuralcategory.png", outrds), 
    width = 11, height = 11, unit = "in", res = 400)
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
p1 <- ggplot(nbrexons %>%
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
p2 <- ggplot(nbrexons, 
             aes(x = structural_category, y = length)) + 
  geom_boxplot() + theme_bw() + 
  facet_wrap(~ dataset, nrow = 1) +
  scale_y_log10() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  xlab("Structural category") + ylab("Inferred transcript length")
cowplot::plot_grid(p0, p1, p2, ncol = 1, labels = c("A", "B", "C"), 
                   rel_heights = c(1.7, 1.08, 0.9), align = "v", axis = "l")
dev.off()

png(gsub("\\.rds$", "_dist_to_tss_tts_by_structuralcategory.png", outrds), 
    width = 11, height = 11, unit = "in", res = 400)
p0 <- ggplot(nbrexons %>% dplyr::filter(!is.na(diff_to_TSS)) %>%
               droplevels(), 
             aes(x = structural_category, y = diff_to_TSS)) + 
  geom_violin() + theme_bw() + 
  facet_wrap(~ dataset, nrow = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  xlab("Structural category") + ylab("Distance to TSS")
p1 <- ggplot(nbrexons %>% dplyr::filter(!is.na(diff_to_TTS)) %>%
               droplevels(), 
             aes(x = structural_category, y = diff_to_TTS)) + 
  geom_violin() + theme_bw() + 
  facet_wrap(~ dataset, nrow = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  xlab("Structural category") + ylab("Distance to TTS")
p2 <- ggplot(nbrexons, 
             aes(x = structural_category, y = perc_A_downstream_TTS)) + 
  geom_violin() + theme_bw() + 
  facet_wrap(~ dataset, nrow = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  xlab("Structural category") + ylab("Percentage A\ndownstream of TTS")
cowplot::plot_grid(p0, p1, p2, ncol = 1, labels = c("A", "B", "C"), 
                   rel_heights = c(1, 1, 1), align = "v", axis = "l")
dev.off()


saveRDS(NULL, file = outrds)
date()
sessionInfo()
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
## Plot overlap of detected transcripts (class code = or c)
df0 <- do.call(dplyr::bind_rows, lapply(setdiff(datasets, "Illumina"), function(ds) {
  read.delim(paste0(
    ds, "/", ifelse(ds == "Illumina", "stringtie_assembly", "flair"),
    "/", ds, flaircond, "/", ds, flaircond, 
    ifelse(ds == "Illumina", "_stringtie_assembly.gffcompare.",
           "_minimap_genome_s_primary_flair_collapse.isoforms.gffcompare."), 
    ds, flaircond, 
    ifelse(ds == "Illumina", "_stringtie_assembly.gtf.tmap", 
           "_minimap_genome_s_primary_flair_collapse.isoforms.gtf.tmap")), 
    header = TRUE, as.is = TRUE) %>%
    dplyr::filter(class_code %in% c("c", "=")) %>%
    dplyr::select(ref_id) %>%
    dplyr::rename(transcript_id = ref_id) %>%
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
## (identified = class code = or c)
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
## different class codes
nbrexons <- do.call(dplyr::bind_rows, lapply(datasets, function(ds) {
  tmp <- as.data.frame(rtracklayer::import(
    paste0(ds, "/", ifelse(ds == "Illumina", "stringtie_assembly", "flair"),
           "/", ds, flaircond, "/", ds, flaircond, 
           ifelse(ds == "Illumina", "_stringtie_assembly.", 
                  "_minimap_genome_s_primary_flair_collapse.isoforms."), 
           "gffcompare.annotated.gtf")))
  dplyr::full_join(
    tmp %>% 
      dplyr::filter(type == "exon") %>%
      dplyr::group_by(transcript_id) %>%
      dplyr::summarize(nbr_exons = length(type),
                       txLength = sum(width)),
    tmp %>%
      dplyr::filter(type == "transcript") %>%
      dplyr::select(transcript_id, class_code),
    by = "transcript_id"
  ) %>% dplyr::mutate(dataset = remapds[ds])
})) %>%
  dplyr::mutate(nbr_exons = replace(nbr_exons, nbr_exons > 2, ">2"))

## ========================================================================== ##
## Plot distribution of class codes
dfcc <- do.call(dplyr::bind_rows, lapply(datasets, function(ds) {
  tmap <- read.delim(paste0(
    ds, "/", ifelse(ds == "Illumina", "stringtie_assembly", "flair"), 
    "/", ds, flaircond, "/", ds, flaircond, 
    ifelse(ds == "Illumina", "_stringtie_assembly.gffcompare.", 
           "_minimap_genome_s_primary_flair_collapse.isoforms.gffcompare."), 
    ds, flaircond, ifelse(ds == "Illumina", "_stringtie_assembly.",
                          "_minimap_genome_s_primary_flair_collapse.isoforms."),
    "gtf.tmap"), 
    header = TRUE, as.is = TRUE)
  as.data.frame(table(tmap$class_code)) %>% 
    dplyr::mutate(dataset = remapds[ds]) %>%
    dplyr::rename(class_code = Var1) %>%
    dplyr::mutate(class_code = as.character(class_code))
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

png(gsub("\\.rds$", "_nbr_exons_read_length_by_classcode.png", outrds), 
    width = 11, height = 11, unit = "in", res = 400)
p0 <- ggplot(dfcc, aes(x = dataset, y = Freq)) + 
  geom_bar(stat = "identity", position = "fill", aes(fill = class_code)) + theme_bw() + 
  scale_fill_manual(values = structure(colfun(15), names = levels(dfcc$class_code)),
                    name = "Overlap type (class code)") + 
  ylab("Fraction of transcripts") + xlab("") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  scale_y_continuous(expand = c(0, 0, 0.05, 0)) + 
  geom_text(data = dfcc %>% dplyr::group_by(dataset) %>%
              dplyr::summarize(nbrTx = sum(Freq)),
            aes(x = dataset, y = 1, label = nbrTx), 
            vjust = 0, nudge_y = 0.02, size = 3)
p1 <- ggplot(nbrexons %>%
               dplyr::group_by(dataset, class_code, nbr_exons) %>%
               dplyr::tally() %>%
               dplyr::ungroup() %>%
               dplyr::mutate(nbr_exons = factor(nbr_exons, 
                                                levels = c(">2", "2", "1"))), 
             aes(x = class_code, y = n, fill = nbr_exons)) + 
  geom_bar(stat = "identity", position = "fill") + theme_bw() + 
  facet_wrap(~ dataset, nrow = 1) +
  xlab("Class code") + ylab("Fraction of transcripts") + 
  scale_fill_manual(values = c("1" = "red", "2" = "blue", ">2" = "grey"),
                    name = "Number of\nexons") + 
  theme(legend.position = "bottom") + 
  scale_y_continuous(expand = c(0, 0, 0.05, 0))
p2 <- ggplot(nbrexons, 
             aes(x = class_code, y = txLength)) + 
  geom_boxplot() + theme_bw() + 
  facet_wrap(~ dataset, nrow = 1) +
  scale_y_log10() + 
  xlab("Class code") + ylab("Inferred transcript length")
cowplot::plot_grid(p0, p1, p2, ncol = 1, labels = c("A", "B", "C"), 
                   rel_heights = c(1.7, 1.08, 0.9), align = "v", axis = "l")
dev.off()


saveRDS(NULL, file = outrds)
date()
sessionInfo()
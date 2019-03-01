args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2, lib.loc = "/home/charlotte/R/x86_64-pc-linux-gnu-library/3.5")
  library(cowplot)
  library(stringr)
})
source("manuscript_results_figures/Rscripts/remap_sample_names.R")

datasets <- strsplit(datasets, ",")[[1]]
datasets <- c(datasets, "Illumina")
names(datasets) <- datasets
conditions <- strsplit(conditions, ",")[[1]]

muted <- c("#DC050C","#E8601C","#7BAFDE","#1965B0","#B17BA6",
           "#882E72","#F1932D","#F6C141","#F7EE55","#4EB265",
           "#CAEDAB","#777777")
colfun <- grDevices::colorRampPalette(muted)

print(tx2gene)
print(datasets)
print(conditions)
print(outrds)

tx2gene <- readRDS(tx2gene)

## Abundance info
abundanceInfo <- lapply(datasets, function(ds) {
  readRDS(paste0(ds, "/output/", ds, "_all_abundances.rds"))
})

## Keep only the desired abundance measures
abundanceInfo <- lapply(abundanceInfo, function(w) {
  w$gene_abundances <- 
    w$gene_abundances[, grep("count__salmon31|count__salmonminimap2_p0.99|count__wubminimap2|count__featurecountsminimap2primary|tpm__salmon$|tpm__StringTie", 
                             colnames(w$gene_abundances))]
  w
})

genelevel <- Reduce(function(...) dplyr::full_join(..., by = "gene"), 
                    lapply(abundanceInfo, function(ab) {
                      ab$gene_abundances %>% tibble::rownames_to_column("gene")
                    })) %>%
  dplyr::left_join(tx2gene %>% dplyr::select(gene, gene_biotype) %>%
                     dplyr::mutate(gene = gsub("\\.[0-9]+$", "", gene)) %>%
                     dplyr::distinct(), by = "gene")


tmp <- genelevel %>% dplyr::select(-gene) %>% dplyr::group_by(gene_biotype) %>% 
  dplyr::summarize_at(vars(dplyr::matches("__count__|__tpm__")), sum) %>%
  dplyr::mutate_at(vars(dplyr::matches("__count__|__tpm__")), funs(./sum(.))) %>%
  tidyr::gather(method, proportion, -gene_biotype) %>%
  tidyr::separate(method, into = c("sample", "aType", "method"), sep = "__") %>%
  dplyr::mutate(method = gsub("featurecounts", "fC", method)) %>% 
  dplyr::mutate(sample = remap[sample]) %>%
  dplyr::left_join(sample_annotation %>% dplyr::select(sample_remap, dataset, condition),
                   by = c("sample" = "sample_remap")) %>%
  dplyr::filter(condition %in% conditions) %>%
  dplyr::group_by(gene_biotype) %>% 
  dplyr::mutate(maxFrac = max(proportion)) %>%
  dplyr::ungroup()

methodcols <- c(salmon = "#882E72", StringTie = "#882E72", 
                salmonminimap2_p0.99 = "#4EB265",
                fCminimap2primary = "#4EB265", salmon31 = "#4EB265",
                wubminimap2 = "#4EB265")

png(gsub("\\.rds$", "_biotype_prop_average.png", outrds), width = 10, 
    height = 15, unit = "in", res = 400)
p1 <- ggplot(tmp %>% 
               dplyr::mutate(dataset = replace(dataset, dataset != "Illumina", 
                                               "Nanopore")) %>% 
               dplyr::group_by(dataset, gene_biotype, method) %>%
               dplyr::summarize(proportion = mean(proportion),
                                maxFrac = max(maxFrac)) %>%
               dplyr::filter(maxFrac >= 0.001) %>%
               dplyr::mutate(method = factor(
                 method, levels = c("fCminimap2primary", 
                                    "salmonminimap2_p0.99", "salmon31", 
                                    "wubminimap2", "salmon", 
                                    "StringTie"))),
             aes(x = method, y = proportion, color = method, fill = method)) + 
  geom_bar(stat = "identity") +
  theme_bw() + facet_wrap(~ gene_biotype, scales = "free_y") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none",
        strip.text = element_text(size = 7)) + 
  xlab("") + ylab("Average abundance proportion") + 
  scale_color_manual(values = methodcols, name = "") + 
  scale_fill_manual(values = methodcols, name = "")

g0 <- 
  ggplot(tmp %>% dplyr::filter(dataset != "Illumina") %>%
           dplyr::mutate(
             gene_biotype = replace(
               gene_biotype, 
               !(gene_biotype %in% c("lincRNA", "Mt_rRNA",
                                     "processed_pseudogene", 
                                     "protein_coding",
                                     "antisense_RNA",
                                     "Mt_tRNA",
                                     "transcribed_processed_pseudogene")),
               "other")),
         aes(x = method, y = proportion, group = gene_biotype, fill = gene_biotype)) + 
  geom_bar(stat = "identity", position = "fill") + 
  facet_wrap(~ sample, ncol = 4) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "bottom") + 
  xlab("") + ylab("Abundance proportion") + 
  scale_fill_manual(values = colfun(8), name = "")

g1 <- 
  ggplot(tmp %>% dplyr::filter(dataset == "Illumina") %>%
           dplyr::mutate(
             gene_biotype = replace(
               gene_biotype, 
               !(gene_biotype %in% c("lincRNA", "Mt_rRNA",
                                     "processed_pseudogene", 
                                     "protein_coding",
                                     "antisense_RNA",
                                     "Mt_tRNA",
                                     "transcribed_processed_pseudogene")),
               "other")),
         aes(x = method, y = proportion, group = gene_biotype, fill = gene_biotype)) + 
  geom_bar(stat = "identity", position = "fill") + 
  facet_wrap(~ sample, ncol = 4) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "bottom") + 
  xlab("") + ylab("Abundance proportion") + 
  scale_fill_manual(values = colfun(8), name = "")

cowplot::plot_grid(
  p1,
  cowplot::plot_grid(g0 + theme(legend.position = "none"), g1, ncol = 1,
                     rel_heights = c(1, 0.45)),
  ncol = 1, labels = c("A", "B"), rel_heights = c(0.8, 1))
dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()

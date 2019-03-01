args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

lowerfun <- function(data, mapping){
  ggplot(data = data, mapping = mapping) +
    geom_point(alpha = 0.3, size = 0.5)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tximport)
  library(ggbeeswarm)
  library(cowplot)
})
source("manuscript_results_figures/Rscripts/remap_sample_names.R")

datasets <- strsplit(datasets, ",")[[1]]
datasets <- setdiff(datasets, "HEK293RNA")
conditions <- strsplit(conditions, ",")[[1]]

print(datasets)
print(conditions)
print(tx2gene)
print(outrds)

tx2gene <- readRDS(tx2gene)

## Read quantifications
salmonfiles_illumina <- list.files("Illumina/salmon31", 
                                   pattern = "quant.sf", 
                                   recursive = TRUE, 
                                   full.names = TRUE)
names(salmonfiles_illumina) <- paste0("Illumina__", 
                                      basename(gsub("/quant.sf", "", 
                                                    salmonfiles_illumina)))
salmonfiles_illumina <- 
  salmonfiles_illumina[names(salmonfiles_illumina) %in%
                         paste0("Illumina__",
                                sample_annotation$sample_orig[sample_annotation$condition %in% 
                                                                conditions])]
txi_illumina <- tximport(salmonfiles_illumina, type = "salmon", txOut = TRUE)
txig_illumina <- summarizeToGene(txi_illumina, tx2gene = tx2gene[, c("tx", "gene")])

## ========================================================================== ##
## Salmon
## ========================================================================== ##
salmonfiles_nanopore <- unlist(lapply(datasets, function(ds) {
  salmon31 <- list.files(file.path(ds, "salmon31"), pattern = "quant.sf",
                         recursive = TRUE, full.names = TRUE)
  names(salmon31) <- paste0("salmon31__", basename(gsub("/quant.sf", "", 
                                                        salmon31)))
  
  salmon31noclip <- list.files(file.path(ds, "salmon31_aligned_noclippedbases"),
                               pattern = "quant.sf", recursive = TRUE, 
                               full.names = TRUE)
  names(salmon31noclip) <- paste0("salmon31noclip__", basename(gsub("/quant.sf", "",
                                                                    salmon31noclip)))
  
  salmonminimap2_p0.8 <- list.files(file.path(ds, "salmonminimap2"), 
                                              pattern = "quant.sf", recursive = TRUE,
                                              full.names = TRUE)
  names(salmonminimap2_p0.8) <- paste0("salmonminimap2_p0.8__",
                                       basename(gsub("/quant.sf", "", 
                                                     salmonminimap2_p0.8)))
  
  salmonminimap2_p0.99 <- list.files(file.path(ds, "salmonminimap2_p0.99"), 
                                               pattern = "quant.sf", recursive = TRUE,
                                               full.names = TRUE)
  names(salmonminimap2_p0.99) <- paste0("salmonminimap2_p0.99__",
                                        basename(gsub("/quant.sf", "", 
                                                      salmonminimap2_p0.99)))
  
  c(salmon31, salmon31noclip, salmonminimap2_p0.8, salmonminimap2_p0.99)
}))

salmonfiles_nanopore <- 
  salmonfiles_nanopore[gsub("salmon31__|salmon31noclip__|salmonminimap2_p0.8__|salmonminimap2_p0.99__", "", names(salmonfiles_nanopore)) %in%
                         sample_annotation$sample_orig[sample_annotation$condition %in% 
                                                         conditions]]

txi_nanopore <- tximport(salmonfiles_nanopore, type = "salmon", txOut = TRUE)
txig_nanopore <- summarizeToGene(txi_nanopore, tx2gene = tx2gene[, c("tx", "gene")])

txi <- as.data.frame(txi_nanopore$counts) %>%
  tibble::rownames_to_column("tx") %>%
  dplyr::full_join(as.data.frame(txi_illumina$abundance) %>% 
                     tibble::rownames_to_column("tx"), by = "tx") %>%
  as.data.frame() %>%
  tibble::column_to_rownames("tx")
txig <- as.data.frame(txig_nanopore$counts) %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::full_join(as.data.frame(txig_illumina$abundance) %>% 
                     tibble::rownames_to_column("gene"), by = "gene") %>%
  as.data.frame() %>%
  tibble::column_to_rownames("gene")

fixCorDf <- function(cordf) {
  cordf %>%
    tibble::rownames_to_column("sample1") %>%
    tidyr::gather(key = "sample2", value = "correlation", -sample1) %>%
    tidyr::separate(sample1, into = c("method1", "sample1"), sep = "__") %>%
    tidyr::separate(sample2, into = c("method2", "sample2"), sep = "__") %>%
    dplyr::mutate(sample1 = remap[sample1],
                  sample2 = remap[sample2]) %>%
    dplyr::mutate(condition1 = stringr::str_extract(sample1, "wt|srpk"),
                  condition2 = stringr::str_extract(sample2, "wt|srpk")) %>%
    dplyr::filter(method1 == "Illumina" & method2 != "Illumina" & 
                    sample1 != sample2 & condition1 == condition2) %>%
    dplyr::mutate(dataset2 = sapply(strsplit(sample2, "_"), 
                                    .subset, 1))
}

spearman_tx <- fixCorDf(as.data.frame(cor(sqrt(txi), method = "spearman")))

spearman_gene <- fixCorDf(as.data.frame(cor(sqrt(txig), method = "spearman")))

pearson_tx <- fixCorDf(as.data.frame(cor(sqrt(txi), method = "pearson")))

pearson_gene <- fixCorDf(as.data.frame(cor(sqrt(txig), method = "pearson")))

png(gsub("\\.rds$", "_salmon_spearman.png", outrds), width = 9, height = 5, 
    unit = "in", res = 400)
p1 <- ggplot(spearman_tx, aes(x = method2, y = correlation, color = dataset2)) + 
  geom_quasirandom(size = 3, alpha = 0.6) + 
  theme_bw() + ggtitle("Transcript") + 
  ylab("Spearman correlation with\nIllumina (sqrt abundance)") + xlab("") + 
  coord_cartesian(ylim = c(0, 1)) + 
  scale_color_manual(values = ds_colors, name = "") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p2 <- ggplot(spearman_gene, aes(x = method2, y = correlation, color = dataset2)) + 
  geom_quasirandom(size = 3, alpha = 0.6) + 
  theme_bw() + ggtitle("Gene") + 
  ylab("Spearman correlation with\nIllumina (sqrt abundance)") + xlab("") + 
  coord_cartesian(ylim = c(0, 1)) + 
  scale_color_manual(values = ds_colors, name = "") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom")
print(cowplot::plot_grid(
  cowplot::plot_grid(
    p1 + theme(legend.position = "none"), 
    p2 + theme(legend.position = "none"),
    labels = c("A", "B"), rel_widths = c(1, 1), nrow = 1
  ),
  cowplot::get_legend(p2),
  rel_heights = c(1, 0.2), ncol = 1)
)
dev.off()

png(gsub("\\.rds$", "_salmon_pearson.png", outrds), width = 9, height = 5, 
    unit = "in", res = 400)
p1 <- ggplot(pearson_tx, aes(x = method2, y = correlation, color = dataset2)) + 
  geom_quasirandom(size = 3, alpha = 0.6) + 
  theme_bw() + ggtitle("Transcript") + 
  ylab("Pearson correlation with\nIllumina (sqrt abundance)") + xlab("") + 
  coord_cartesian(ylim = c(0, 1)) + 
  scale_color_manual(values = ds_colors, name = "") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p2 <- ggplot(pearson_gene, aes(x = method2, y = correlation, color = dataset2)) + 
  geom_quasirandom(size = 3, alpha = 0.6) + 
  theme_bw() + ggtitle("Gene") + 
  ylab("Pearson correlation with\nIllumina (sqrt abundance)") + xlab("") + 
  coord_cartesian(ylim = c(0, 1)) + 
  scale_color_manual(values = ds_colors, name = "") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom")
print(cowplot::plot_grid(
  cowplot::plot_grid(
    p1 + theme(legend.position = "none"), 
    p2 + theme(legend.position = "none"),
    labels = c("A", "B"), rel_widths = c(1, 1), nrow = 1
  ),
  cowplot::get_legend(p2),
  rel_heights = c(1, 0.2), ncol = 1)
)
dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(ggbeeswarm)
})
source("manuscript_results_figures/Rscripts/remap_sample_names.R")

datasets <- strsplit(datasets, ",")[[1]]
datasets <- c(datasets, "Illumina")
names(datasets) <- datasets
conditions <- strsplit(conditions, ",")[[1]]

print(tx2gene)
print(datasets)
print(conditions)
print(outrds)

## Define base colors
muted <- c("#DC050C","#E8601C","#7BAFDE","#1965B0","#B17BA6",
           "#882E72","#F1932D","#F6C141","#F7EE55","#4EB265",
           "#90C987","#CAEDAB","#777777")
colfun <- grDevices::colorRampPalette(muted)

abundances <- lapply(datasets, function(ds) {
  readRDS(paste0(ds, "/output/", ds, "_all_abundances.rds"))
})

## Keep only the desired abundance measures
abundances <- lapply(abundances, function(w) {
  w$tx_abundances <- 
    w$tx_abundances[, grep("count__salmon31|count__salmonminimap2_p0.99|count__wubminimap2|tpm__salmon$|tpm__StringTie", 
                           colnames(w$tx_abundances))]
  w$gene_abundances <- 
    w$gene_abundances[, grep("count__salmon31|count__salmonminimap2_p0.99|count__wubminimap2|count__featurecountsminimap2primary|tpm__salmon$|tpm__StringTie", 
                             colnames(w$gene_abundances))]
  w
})

tx2gene <- readRDS(tx2gene)
## Remove version information from tx2gene
tx2gene <- tx2gene %>% dplyr::mutate(tx = gsub("\\.[0-9]+$", "", tx),
                                     gene = gsub("\\.[0-9]+$", "", gene))

## Transcript-level
txlevel <- Reduce(function(...) dplyr::full_join(..., by = "tx"), 
                  lapply(abundances, function(ab) {
                    ab$tx_abundances %>% tibble::rownames_to_column("tx")
                  })) %>%
  tibble::column_to_rownames("tx")

## Gene-level
genelevel <- Reduce(function(...) dplyr::full_join(..., by = "gene"), 
                  lapply(abundances, function(ab) {
                    ab$gene_abundances %>% tibble::rownames_to_column("gene")
                  })) %>%
  tibble::column_to_rownames("gene")

cortx <- list()
corgene <- list()

## Help functions
calcCor <- function(abundancemat, method) {
  cor(abundancemat,
      method = method, use = "pairwise.complete.obs") %>%
    as.data.frame() %>% tibble::rownames_to_column("sample1") %>%
    tidyr::gather(key = "sample2", value = "correlation", -sample1) %>%
    tidyr::separate(sample1, into = c("sample1", "type1", "method1"), sep = "__") %>%
    tidyr::separate(sample2, into = c("sample2", "type2", "method2"), sep = "__") %>%
    tidyr::unite(col = method1, type1, method1, sep = "__") %>%
    tidyr::unite(col = method2, type2, method2, sep = "__") %>%
    dplyr::mutate(sample1 = remap[sample1], 
                  sample2 = remap[sample2]) %>%
    dplyr::mutate(dataset1 = sapply(strsplit(sample1, "_"), .subset, 1),
                  dataset2 = sapply(strsplit(sample2, "_"), .subset, 1),
                  condition1 = sapply(strsplit(sample1, "_"), .subset, 2),
                  condition2 = sapply(strsplit(sample2, "_"), .subset, 2)) %>% 
    dplyr::filter(condition1 %in% conditions & condition2 %in% conditions) %>%
    dplyr::mutate(dtype = ifelse(
      dataset1 == "Illumina" & dataset2 == "Illumina", "Illumina", 
      ifelse((dataset1 == "Illumina" & dataset2 != "Illumina") | 
               (dataset1 != "Illumina" & dataset2 == "Illumina"),
             "Illumina-Nanopore", "Nanopore"))) %>%
    dplyr::mutate(combination = ifelse(dataset1 == dataset2, "same dataset", 
                                       "different datasets"))
}

plotCor <- function(cordf, ylab, title) {
  tmp <- cordf %>% 
    dplyr::filter(dataset1 != "ONT-RNA001-HEK" & dataset2 != "ONT-RNA001-HEK") %>%
    dplyr::mutate(method1 = gsub("featurecounts", "fC", 
                                 gsub("count__", "", method1)),
                  method2 = gsub("featurecounts", "fC", 
                                 gsub("count__", "", method2))) %>%
    dplyr::filter(gsub("tpm__", "", gsub("31", "", method1)) == 
                    gsub("tpm__", "", gsub("31", "", method2)) & 
                    sample1 != sample2 & 
                    condition1 == condition2 & 
                    (dataset1 < dataset2 | (dataset1 == dataset2 & 
                                              sample1 < sample2))) %>%
    dplyr::mutate(
      method1 = replace(
        method1, dtype == "Illumina-Nanopore", 
        gsub("tpm__", "", gsub("31", "", 
                               method1[dtype == "Illumina-Nanopore"]))),
      method2 = replace(
        method2, dtype == "Illumina-Nanopore", 
        gsub("tpm__", "", gsub("31", "", 
                               method2[dtype == "Illumina-Nanopore"])))
    )
  ncol <- length(unique(interaction(tmp$dataset1, tmp$dataset2)))
  ggplot(tmp, 
         aes(x = method1, y = correlation, color = interaction(dataset1, dataset2),
             shape = combination)) + 
    geom_quasirandom(size = 3, alpha = 0.6) + theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1, size = 12)) + 
    xlab("") + ylab(ylab) + 
    scale_color_manual(values = colfun(ncol), name = "") + 
    expand_limits(y = c(0, 1)) + scale_shape_discrete(name = "") + 
    facet_grid(~ dtype, scales = "free_x", space = "free_x") + 
    ggtitle(title)
}

plotCorBetween <- function(cordf, ylab, title) {
  ggplot(cordf %>% 
           dplyr::mutate(method1 = gsub("featurecounts", "fC", 
                                        gsub("count__", "", method1)),
                         method2 = gsub("featurecounts", "fC", 
                                        gsub("count__", "", method2))) %>%
           dplyr::filter(method1 < method2 & sample1 == sample2) %>%
           dplyr::mutate(dataset1 = factor(dataset1, levels = 
                                             ds_order[ds_order %in% dataset1])), 
         aes(x = interaction(method1, method2), y = correlation, color = dataset1)) + 
    geom_quasirandom(size = 4, alpha = 0.8) + theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
    xlab("") + ylab(ylab) + ggtitle(title) + 
    scale_color_manual(values = ds_colors, 
                       name = "") + expand_limits(y = c(0, 1)) + 
    facet_grid(~ dtype, scales = "free_x", space = "free_x")
}

for (m in c("Pearson", "Spearman")) {
  cortx[[m]] <- calcCor(sqrt(txlevel), method = tolower(m))
  corgene[[m]] <- calcCor(sqrt(genelevel), method = tolower(m))

  ptx <- plotCor(cordf = cortx[[m]], 
                 ylab = paste0(m, " correlation between\nreplicate pairs, ", 
                               "sqrt(abundances)"),
                 title = "Transcript")
  pg <- plotCor(cordf = corgene[[m]],
                ylab = paste0(m, " correlation between\nreplicate pairs, ", 
                              "sqrt(abundances)"),
                title = "Gene")
  
  png(gsub("\\.rds", paste0("_", m, ".png"), outrds), 
      width = 18, height = 6, unit = "in", res = 400)
  print(cowplot::plot_grid(
    ptx + theme(legend.position = "none"),
    pg + theme(legend.position = "none"),
    cowplot::get_legend(pg),
    nrow = 1, rel_widths = c(1, 1.15, 0.5), align = "h", axis = "t",
    labels = c("A", "B", "")
  ))
  dev.off()
  
  ptxb <- plotCorBetween(cordf = cortx[[m]],
                         ylab = paste0(m, " correlation between\nmethod pairs, ", 
                                       "sqrt(abundances)"),
                         title = "Transcript")
  pgb <- plotCorBetween(cordf = corgene[[m]],
                        ylab = paste0(m, " correlation between\nmethod pairs, ", 
                                      "sqrt(abundances)"),
                        title = "Gene")
  
  png(gsub("\\.rds", paste0("_betweenmethods_", m, ".png"), outrds), 
      width = 13, height = 6, unit = "in", res = 400)
  print(cowplot::plot_grid(
    cowplot::plot_grid(
      ptxb + theme(legend.position = "none"),
      pgb + theme(legend.position = "none"),
      nrow = 1, rel_widths = c(1, 1), labels = c("A", "B"),
      align = "h", axis = "tb"),
    cowplot::get_legend(pgb), nrow = 1, rel_widths = c(2, 0.4),
    labels = c("", ""), align = "h", axis = "t"
  ))
  dev.off()
}

saveRDS(NULL, file = outrds)
date()
sessionInfo()



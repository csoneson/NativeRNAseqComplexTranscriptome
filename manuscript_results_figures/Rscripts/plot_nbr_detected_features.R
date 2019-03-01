args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(stringr)
})
source("manuscript_results_figures/Rscripts/remap_sample_names.R")

datasets <- strsplit(datasets, ",")[[1]]
datasets <- c(datasets, "Illumina")
names(datasets) <- datasets
conditions <- strsplit(conditions, ",")[[1]]

print(datasets)
print(conditions)
print(tx2gene)
print(outrds)

tx2gene <- readRDS(tx2gene)

## Abundance info
abundanceInfo <- lapply(datasets, function(ds) {
  readRDS(paste0(ds, "/output/", ds, "_all_abundances.rds"))
})

## Keep only the desired abundance measures
abundanceInfo <- lapply(abundanceInfo, function(w) {
  w$tx_abundances <- 
    w$tx_abundances[, grep("count__salmon31|count__salmonminimap2_p0.99|count__wubminimap2|count__salmon$|count__StringTie", 
                           colnames(w$tx_abundances))]
  w$gene_abundances <- 
    w$gene_abundances[, grep("count__salmon31|count__salmonminimap2_p0.99|count__wubminimap2|count__featurecountsminimap2primary|count__salmon$|count__StringTie", 
                             colnames(w$gene_abundances))]
  w
})

## Abundance matrices
txlevel <- Reduce(function(...) dplyr::full_join(..., by = "tx"), 
                  lapply(abundanceInfo, function(ab) {
                    ab$tx_abundances %>% tibble::rownames_to_column("tx")
                  })) %>%
  tibble::column_to_rownames("tx")
txlevel[is.na(txlevel)] <- 0
txlevel <- round(txlevel)
genelevel <- Reduce(function(...) dplyr::full_join(..., by = "gene"), 
                    lapply(abundanceInfo, function(ab) {
                      ab$gene_abundances %>% tibble::rownames_to_column("gene")
                    })) %>%
  tibble::column_to_rownames("gene")
genelevel[is.na(genelevel)] <- 0
genelevel <- round(genelevel)

png(gsub("\\.rds$", "_cumulative_abundance_genes.png", outrds), width = 7, 
    height = 7, unit = "in", res = 400)
gl <- data.frame(apply(genelevel[rowSums(genelevel) > 0, ], 
                       2, function(w) cumsum(sort(w, decreasing = TRUE))/sum(w)), 
                 check.names = FALSE) %>%
  dplyr::mutate(x = seq_len(sum(rowSums(genelevel) > 0))) %>%
  dplyr::mutate(x = x/max(x)) %>% 
  tidyr::gather("method", "cumulativeAbundance", -x) %>%
  tidyr::separate(method, into = c("sample", "aType", "method"), sep = "__") %>%
  dplyr::mutate(sample = remap[sample]) %>%
  dplyr::left_join(sample_annotation %>% dplyr::select(sample_remap, condition, dataset), 
                   by = c("sample" = "sample_remap")) %>%
  dplyr::filter(condition %in% conditions) %>%
  dplyr::filter(method %in% c("salmon31", "salmon"))
ggplot(gl, aes(x = x, y = cumulativeAbundance, group = sample, color = dataset)) + 
  geom_line(size = 1.1) + coord_cartesian(xlim = c(0, 0.25)) + 
  theme_bw() + xlab("Fraction of top expressed genes") + 
  ylab("Cumulative abundance") + 
  theme(legend.position = "bottom") + 
  scale_color_manual(values = ds_colors, name = "") + 
  guides(color = guide_legend(override.aes = list(size = 2))) 
dev.off()

## Number of detected genes/transcripts and number of reads
totGenes <- 
  rbind(
    data.frame(totFeatures = colSums(genelevel >= 1), 
               method = paste0(colnames(genelevel), "__totFeatures"), 
               stringsAsFactors = FALSE),
    data.frame(totFeatures = colSums(genelevel),
               method = paste0(colnames(genelevel), "__totReads"),
               stringsAsFactors = FALSE)
  ) %>%
  tidyr::separate(method, into = c("sample", "type", "method", "fType"), sep = "__") %>%
  dplyr::mutate(sample = remap[sample])%>% 
  dplyr::left_join(sample_annotation %>% dplyr::select(sample_remap, dataset, condition),
                   by = c("sample" = "sample_remap")) %>%
  dplyr::filter(condition %in% conditions) %>%
  dplyr::select(sample, dataset, method, fType, totFeatures) %>%
  dplyr::mutate(
    method = replace(method, method == "featurecountsminimap2primary",
                     "Number of detected genes, fCminimap2primary"),
    method = replace(method, method == "salmon31",
                     "Number of detected genes, salmon31"),
    method = replace(method, method == "salmonminimap2_p0.99",
                     "Number of detected genes, salmonminimap2_p0.99"),
    method = replace(method, method == "wubminimap2",
                     "Number of detected genes, wubminimap2"),
    method = replace(method, method == "salmon",
                     "Number of detected genes, salmon31"),
    method = replace(method, method == "StringTie",
                     "Number of detected genes, StringTie")
  ) %>%
  dplyr::mutate(
    method = factor(method, 
                    levels = c("Number of detected genes, salmonminimap2_p0.99",
                               "Number of detected genes, fCminimap2primary",
                               "Number of detected genes, salmon31",
                               "Number of detected genes, wubminimap2",
                               "Number of detected genes, StringTie")
    )
  )

totTx <- 
  rbind(
    data.frame(totFeatures = colSums(txlevel >= 1), 
               method = paste0(colnames(txlevel), "__totFeatures"), 
               stringsAsFactors = FALSE),
    data.frame(totFeatures = colSums(txlevel),
               method = paste0(colnames(txlevel), "__totReads"),
               stringsAsFactors = FALSE)
  ) %>%
  tidyr::separate(method, into = c("sample", "type", "method", "fType"), sep = "__") %>%
  dplyr::mutate(sample = remap[sample])%>% 
  dplyr::left_join(sample_annotation %>% dplyr::select(sample_remap, dataset, condition),
                   by = c("sample" = "sample_remap")) %>%
  dplyr::filter(condition %in% conditions) %>%
  dplyr::select(sample, dataset, method, fType, totFeatures) %>%
  dplyr::mutate(
    method = replace(method, method == "salmon31",
                     "Number of detected transcripts, salmon31"),
    method = replace(method, method == "salmonminimap2_p0.99",
                     "Number of detected transcripts, salmonminimap2_p0.99"),
    method = replace(method, method == "wubminimap2",
                     "Number of detected transcripts, wubminimap2"),
    method = replace(method, method == "salmon",
                     "Number of detected transcripts, salmon31"),
    method = replace(method, method == "StringTie",
                     "Number of detected transcripts, StringTie")
  ) %>%
  dplyr::mutate(
    method = factor(method, 
                    levels = c("Number of detected transcripts, salmonminimap2_p0.99",
                               "Number of detected transcripts, salmon31",
                               "Number of detected transcripts, wubminimap2",
                               "Number of detected transcripts, StringTie")
    )
  )

totals <- list(transcripts = totTx, 
               genes = totGenes)

methodcols <- c(`salmonminimap2_p0.99` = "#1965B0",
                `fCminimap2primary` = "#882E72",
                `salmon31` = "#DC050C",
                `wubminimap2` = "#55A1B1",
                `StringTie` = "darkgrey")

for (tl in names(totals)) {
  png(gsub("\\.rds$", paste0("_nbr_detected_", tl, "_vs_nbr_reads.png"), outrds),
      width = 7, height = 8, unit = "in", res = 400) 
  print(ggplot(totals[[tl]] %>% dplyr::filter(dataset != "Illumina") %>%
                 tidyr::spread(key = fType, value = totFeatures),
               aes(x = totReads, y = totFeatures, group = sample, 
                   color = gsub(paste0("Number of detected ", tl, ", "), "", method))) + 
          geom_line(color = "grey") + geom_point(size = 2) + 
          facet_wrap(~ dataset, scales = "free_x") + theme_bw() + 
          scale_color_manual(values = methodcols, name = "") + 
          theme(legend.position = "bottom") + 
          xlab("Number of assigned reads") + 
          ylab(paste0("Number of detected ", tl)))
  dev.off()
}

gglayers <- list(
  geom_bar(stat = "identity", position = "dodge"),
  theme_bw(),
  scale_y_continuous(expand = c(0, 0, 0.05, 0)),
  facet_grid(~ dataset, scales = "free_x", space = "free_x"),
  scale_fill_manual(values = methodcols, name = ""),
  theme(legend.position = "bottom",
        strip.text = element_text(size = 8)),
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)),
  xlab("")
)
pt <- ggplot(totals[["transcripts"]] %>% dplyr::filter(fType == "totFeatures") %>%
               dplyr::mutate(sample = removeDatasetFromSample(sample, dataset)),
             aes(x = sample, y = totFeatures, 
                 fill = gsub(paste0("Number of detected transcripts, "), "", method))) + 
  gglayers + 
  ylab(paste0("Number of detected features")) + ggtitle("Transcripts")
pg <- ggplot(totals[["genes"]] %>% dplyr::filter(fType == "totFeatures") %>%
               dplyr::mutate(sample = removeDatasetFromSample(sample, dataset)),
             aes(x = sample, y = totFeatures, 
                 fill = gsub(paste0("Number of detected genes, "), "", method))) +
  gglayers + 
  ylab(paste0("Number of detected features")) + ggtitle("Genes")

################################################################################
## Subsampling - saturation
################################################################################
txCountMatrix <- as.matrix(txlevel)
geneCountMatrix <- as.matrix(genelevel)

## Remap column names and add a column for each data set/condition with the sum
## of the counts across all the samples therein
remap_and_add_sum <- function(countmatrix, sample_annotation, conditions) {
  tmp <- data.frame(countmatrix, check.names = FALSE) %>%
    tibble::rownames_to_column("feature") %>%
    tidyr::gather(key = "sample", value = "count", -feature) %>%
    tidyr::separate(sample, into = c("sample", "mtype", "method"), sep = "__") %>%
    dplyr::mutate(sample = remap[sample]) %>%
    dplyr::left_join(sample_annotation %>% dplyr::select(sample_remap, dataset, condition),
                     by = c("sample" = "sample_remap")) %>%
    dplyr::filter(condition %in% conditions)
  tmpds <- tmp %>% dplyr::group_by(feature, mtype, method, dataset, condition) %>%
    dplyr::summarize(count = sum(count)) %>% 
    dplyr::mutate(sample = paste0(dataset, "_", condition, "_sum"))
  dplyr::bind_rows(tmp, tmpds) %>%
    dplyr::select(-dataset, -condition) %>%
    tidyr::unite(col = "sample", sample, mtype, method, sep = "__") %>%
    tidyr::spread(sample, count) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("feature") %>%
    as.matrix()
}
txCountMatrix <- remap_and_add_sum(txCountMatrix, sample_annotation, conditions)
geneCountMatrix <- remap_and_add_sum(geneCountMatrix, sample_annotation, conditions)

subfracs <- c(seq(from = 0.001, to = 0.025, length.out = 10), 
              c(0.035, 0.045, 0.06, 0.08), 
              seq(from = 0.1, to = 1, by = 0.05))

subsampling <- do.call(rbind, lapply(subfracs, function(i) {
  print(i)
  txmat <- DropletUtils::downsampleMatrix(txCountMatrix, prop = i, bycol = TRUE)
  genemat <- DropletUtils::downsampleMatrix(geneCountMatrix, prop = i, bycol = TRUE)
  dplyr::full_join(
    data.frame(sample = colnames(txCountMatrix), 
               nTx = colSums(txmat >= 1),
               nReadsTx = colSums(txmat),
               stringsAsFactors = FALSE),
    data.frame(sample = colnames(geneCountMatrix),
               nGene = colSums(genemat >= 1),
               nReadsGene = colSums(genemat),
               stringsAsFactors = FALSE)
  ) %>% tidyr::separate(sample, into = c("sample", "cType", "method"), sep = "__") %>%
    tidyr::separate(sample, into = c("dataset", "condition", "replicate"), 
                    sep = "_", remove = FALSE) %>%
    dplyr::mutate(sampleFrac = i) %>%
    dplyr::select(-condition)
})) 

png(gsub("\\.rds$", "_transcript_detection_saturation.png", outrds),
    width = 15, height = 6, unit = "in", res = 400)
q2 <- ggplot(subsampling %>% dplyr::filter(nReadsTx < 1.1e6) %>%
               dplyr::filter(replicate != "sum") %>% 
               dplyr::filter(dataset == "Illumina" | method == "salmonminimap2_p0.99"), 
             aes(x = nReadsTx, y = nTx, group = sample, color = dataset)) + 
  geom_line(color = "black") + geom_point() + 
  scale_color_manual(values = ds_colors, name = "") + 
  theme_bw() + xlab("Sampled number of assigned reads") + 
  theme(legend.position = "bottom") + ggtitle("Transcripts") + 
  ylab("Number of detected features")
q3 <- ggplot(subsampling %>% dplyr::filter(nReadsGene < 1.1e6) %>%
               dplyr::filter(replicate != "sum") %>% 
               dplyr::filter(dataset == "Illumina" | method == "salmonminimap2_p0.99"), 
             aes(x = nReadsGene, y = nGene, group = sample, color = dataset)) + 
  geom_line(color = "black") + geom_point() + 
  scale_color_manual(values = ds_colors, name = "") + 
  theme_bw() + xlab("Sampled number of assigned reads") + 
  theme(legend.position = "bottom") + ggtitle("Genes") + 
  ylab("Number of detected features")
q23 <- cowplot::plot_grid(q2 + theme(legend.position = "none"),
                          q3 + theme(legend.position = "none"),
                          nrow = 1, labels = c("D", "E"))
q23l <- cowplot::plot_grid(
  q23, cowplot::get_legend(q2 + guides(colour = guide_legend(override.aes = 
                                                               list(size = 4)))), 
  ncol = 1, rel_heights = c(1, 0.15), labels = "")
print(q23l)
dev.off()

dsrep_colors <- structure(rep("black", length(unique(paste0(subsampling$dataset, subsampling$replicate)))), names = unique(paste0(subsampling$dataset, subsampling$replicate)))
idx <- grep("sum$", names(dsrep_colors), value = TRUE)
dsrep_colors[idx] <- ds_colors[gsub("sum$", "", idx)]
png(gsub("\\.rds$", "_transcript_detection_saturation_byds.png", outrds),
    width = 15, height = 6, unit = "in", res = 400)
q2ds <- ggplot(subsampling %>%
                 dplyr::filter(dataset == "Illumina" | 
                                 method == "salmonminimap2_p0.99") %>% 
                 dplyr::group_by(dataset) %>%
                 dplyr::filter(nReadsTx < max(nReadsTx[replicate != "sum"])) %>%
                 dplyr::ungroup(), 
             aes(x = nReadsTx, y = nTx, group = sample)) + 
  geom_line(aes(color = paste0(dataset, replicate), size = (replicate == "sum"))) + 
  facet_wrap(~ dataset, scales = "free", nrow = 1) + 
  scale_color_manual(values = dsrep_colors, name = "") + 
  scale_size_manual(values = c(`TRUE` = 1.25, `FALSE` = 1)) + 
  theme_bw() + xlab("Sampled number of assigned reads") + 
  theme(legend.position = "bottom") + ggtitle("Transcripts") + 
  ylab("Number of detected features")
q3ds <- ggplot(subsampling %>%
                 dplyr::filter(dataset == "Illumina" | 
                                 method == "salmonminimap2_p0.99") %>% 
                 dplyr::group_by(dataset) %>%
                 dplyr::filter(nReadsGene < max(nReadsGene[replicate != "sum"])) %>%
                 dplyr::ungroup(), 
               aes(x = nReadsGene, y = nGene, group = sample)) + 
  geom_line(aes(color = paste0(dataset, replicate), size = (replicate == "sum"))) + 
  facet_wrap(~ dataset, scales = "free", nrow = 1) + 
  scale_color_manual(values = dsrep_colors, name = "") + 
  scale_size_manual(values = c(`TRUE` = 1.25, `FALSE` = 1)) + 
  theme_bw() + xlab("Sampled number of assigned reads") + 
  theme(legend.position = "bottom") + ggtitle("Genes") + 
  ylab("Number of detected features")
q23ds <- cowplot::plot_grid(q2ds + theme(legend.position = "none"),
                            q3ds + theme(legend.position = "none"),
                            ncol = 1, labels = c("A", "B"))
print(q23ds)
dev.off()

################################################################################
## Detection rate by transcript length
txl <- txlevel %>% 
  tibble::rownames_to_column("tx") %>%
  tidyr::gather(key = sample, value = abundance, -tx) %>%
  tidyr::separate(sample, into = c("sample", "dType", "method"), sep = "__") %>%
  dplyr::mutate(sample = remap[sample]) %>%
  dplyr::left_join(sample_annotation %>% dplyr::select(sample_remap, dataset, 
                                                       condition),
                   by = c("sample" = "sample_remap")) %>%
  dplyr::filter(condition %in% conditions) %>%
  dplyr::left_join(tx2gene %>% dplyr::select(tx, txlength) %>%
                     dplyr::mutate(tx = gsub("\\.[0-9]+$", "", tx)),
                   by = "tx") %>%
  dplyr::mutate(lengthbin = Hmisc::cut2(txlength, 
                cuts = c(0, 400, 800, 1200, 1600, 2000, 2400, 2800, 
                         3200, 3600, max(txlength + 1)))) %>%
  dplyr::filter(method == "salmonminimap2_p0.99" | (dataset == "Illumina" & 
                                                      method == "salmon")) %>%
  dplyr::group_by(sample, dataset, lengthbin) %>%
  dplyr::summarize(n = length(abundance),
                   detRate = mean(abundance >= 1)) %>%
  dplyr::group_by(dataset, lengthbin) %>%
  dplyr::summarize(detRate = mean(detRate), 
                   n = mean(n))

png(gsub("\\.rds$", "_detrate_by_length.png", outrds), 
    width = 12, height = 4, unit = "in", res = 400)
drbl <- ggplot(txl, aes(x = lengthbin, y = detRate)) + 
  geom_bar(aes(fill = dataset), stat = "identity", position = "dodge") + theme_bw() + 
  scale_fill_manual(values = ds_colors, name = "") + 
  facet_wrap(~ dataset, nrow = 1) + 
  xlab("Transcript length") + ylab("Transcript detection rate") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none")
print(drbl)
dev.off()

################################################################################

png(gsub("\\.rds$", paste0("_nbr_detected_features.png"), outrds), 
    width = 14, height = 18, unit = "in", res = 400)
print(
  cowplot::plot_grid(
    cowplot::plot_grid(
      pt + theme(legend.position = "none"), 
      pg,
      drbl,
      labels = c("A", "B", "C"), rel_heights = c(1, 1.2, 1.15), ncol = 1,
      align = "v", axis = "lr"
    ),
    q23l, labels = c("", ""), rel_heights = c(4, 1.75), ncol = 1
  )
)
dev.off()


saveRDS(NULL, file = outrds)
date()
sessionInfo()

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
})

print(abundances)
print(outrds)

abundances <- readRDS(abundances)

## ========================================================================== ##
## Transcript-level
txab <- abundances$tx_abundances

## Get quant methods (only for counts)
quantmethods <- unique(sapply(strsplit(colnames(txab), "__"), .subset, 3)[sapply(strsplit(colnames(txab), "__"), .subset, 2) == "count"])

edger_dge_tx <- lapply(quantmethods, function(qm) {
  pdf(paste0(dirname(outrds), "/dge_edger_plots/", 
             gsub("\\.rds$", paste0("_transcript_", qm, ".pdf"), basename(outrds))))
  
  counts <- txab[, grep(paste0("__count__", qm, "$"), colnames(txab))]
  grp <- c(str_match(colnames(counts), "srpk|wt|WT|Srpk|Sprk|sprk"))
  meta <- data.frame(sample = colnames(counts), 
                     group = grp,
                     stringsAsFactors = FALSE)
  print(meta)
  
  dge <- DGEList(counts = counts, samples = meta, group = meta$group)
  dge <- calcNormFactors(dge)
  design <- model.matrix(~ group, data = dge$samples)
  dge <- estimateDisp(dge, design = design)
  mds <- plotMDS(dge, top = 100, ndim = 2, gene.selection = "common", 
                 plot = FALSE)$cmdscale.out
  colnames(mds) <- paste0("MDS", seq_len(ncol(mds)))
  mds <- as.data.frame(mds) %>% tibble::rownames_to_column(var = "sample") %>%
    dplyr::full_join(dge$samples) %>% 
    dplyr::mutate(sample = gsub(paste0("__", qm), "", sample))
  print(ggplot(mds, aes(x = MDS1, y = MDS2, label = sample, color = group)) + 
          geom_point(size = 5) + geom_label_repel(size = 2) +
          ggtitle(paste0(qm, ", transcript")))
  plotBCV(dge, main = paste0(qm, ", transcript"))
  fit <- glmQLFit(dge, design = design)
  qlf <- glmQLFTest(fit, coef = 2)
  tt <- topTags(qlf, n = Inf, sort.by = "none")
  plotSmear(qlf, de.tags = rownames(tt)[tt$FDR <= 0.05], 
            main = paste0(qm, ", transcript"))
  dev.off()
  tt$table %>% tibble::rownames_to_column("tx") %>%
    dplyr::mutate(quantmethod = qm)
})

## ========================================================================== ##
## Gene-level
geneab <- abundances$gene_abundances

## Get quant methods (only for counts)
quantmethods <- unique(sapply(strsplit(colnames(geneab), "__"), .subset, 3)[sapply(strsplit(colnames(geneab), "__"), .subset, 2) == "count"])

edger_dge_gene <- lapply(quantmethods, function(qm) {
  pdf(paste0(dirname(outrds), "/dge_edger_plots/", 
             gsub("\\.rds$", paste0("_gene_", qm, ".pdf"), basename(outrds))))
  
  counts <- geneab[, grep(paste0("__count__", qm, "$"), colnames(geneab))]
  grp <- c(str_match(colnames(counts), "srpk|wt|WT|Srpk|Sprk|sprk"))
  meta <- data.frame(sample = colnames(counts), 
                     group = grp,
                     stringsAsFactors = FALSE)
  print(meta)
  
  dge <- DGEList(counts = counts, samples = meta, group = meta$group)
  dge <- calcNormFactors(dge)
  design <- model.matrix(~ group, data = dge$samples)
  dge <- estimateDisp(dge, design = design)
  mds <- plotMDS(dge, top = 100, ndim = 2, gene.selection = "common", 
                 plot = FALSE)$cmdscale.out
  colnames(mds) <- paste0("MDS", seq_len(ncol(mds)))
  mds <- as.data.frame(mds) %>% tibble::rownames_to_column(var = "sample") %>%
    dplyr::full_join(dge$samples) %>% 
    dplyr::mutate(sample = gsub(paste0("__", qm), "", sample))
  print(ggplot(mds, aes(x = MDS1, y = MDS2, label = sample, color = group)) + 
          geom_point(size = 5) + geom_label_repel(size = 2) +
          ggtitle(paste0(qm, ", gene")))
  plotBCV(dge, main = paste0(qm, ", gene"))
  fit <- glmQLFit(dge, design = design)
  qlf <- glmQLFTest(fit, coef = 2)
  tt <- topTags(qlf, n = Inf, sort.by = "none")
  plotSmear(qlf, de.tags = rownames(tt)[tt$FDR <= 0.05], 
            main = paste0(qm, ", gene"))
  dev.off()
  tt$table %>% tibble::rownames_to_column("gene") %>%
    dplyr::mutate(quantmethod = qm)
})

saveRDS(list(edger_dge_tx = edger_dge_tx, edger_dge_gene = edger_dge_gene),
        file = outrds)

date()
sessionInfo()

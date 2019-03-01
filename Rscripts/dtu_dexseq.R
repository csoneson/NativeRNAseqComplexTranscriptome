args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(DEXSeq)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
})

print(abundances)
print(tx2gene)
print(outrds)

abundances <- readRDS(abundances)
tx2gene <- readRDS(tx2gene)
tx2gene <- tx2gene %>% dplyr::mutate(tx = gsub("\\.[0-9]+$", "", tx),
                                     gene = gsub("\\.[0-9]+$", "", gene))

txab <- abundances$tx_abundances

## Get quant methods (only for counts)
quantmethods <- unique(sapply(strsplit(colnames(txab), "__"), .subset, 3)[sapply(strsplit(colnames(txab), "__"), .subset, 2) == "count"])
names(quantmethods) <- quantmethods

dexseq_dtu <- lapply(quantmethods, function(qm) {
  counts <- txab[, grep(paste0("__count__", qm, "$"), colnames(txab))]
  grp <- c(str_match(colnames(counts), "srpk|wt|WT|Srpk|Sprk|sprk"))
  meta <- data.frame(sample = colnames(counts), 
                     condition = factor(grp),
                     stringsAsFactors = FALSE)
  print(meta)
  
  counts <- counts[rowSums(counts) > ncol(counts)/2, ]

  BPPARAM <- MulticoreParam(workers = 8)
  
  ## Generate DEXSeqDataSet
  message("Generating DEXSeqDataSet...")
  dxd <- DEXSeqDataSet(countData = round(counts), 
                       sampleData = meta,
                       design = ~ sample + exon + condition:exon, 
                       featureID = rownames(counts),
                       groupID = tx2gene$gene[match(rownames(counts), tx2gene$tx)])
  
  ## Estimate size factors and dispersions
  dxd <- estimateSizeFactors(dxd)
  message("Estimating dispersions...")
  dxd <- estimateDispersions(dxd, BPPARAM = BPPARAM)
  
  ## Run the test and estimate fold changes
  message("Running test...")
  dxd <- testForDEU(dxd, BPPARAM = BPPARAM)
  dxd <- estimateExonFoldChanges(dxd, BPPARAM = BPPARAM)
  print(dim(dxd))
  
  message("Summarizing results on gene level...")
  res <- DEXSeqResults(dxd)
  pgq <- perGeneQValue(res, p = "pvalue")
  res <- as.data.frame(res)
  res$genomicData <- NULL
  res$padj_gene <- pgq[match(res$groupID, names(pgq))]
  
  res %>% dplyr::rename(gene = groupID, tx = featureID) %>%
    dplyr::left_join(tx2gene %>% dplyr::select(tx, symbol, gene_biotype) %>%
                       dplyr::distinct()) %>%
    dplyr::select(gene, tx, symbol, gene_biotype, everything()) %>%
    dplyr::mutate(quantmethod = qm)
})

saveRDS(dexseq_dtu, file = outrds)

date()
sessionInfo()

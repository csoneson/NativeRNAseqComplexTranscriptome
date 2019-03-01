args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(salmondir)
print(metafile)
print(tx2gene)
print(outrds)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(DEXSeq))

## List Salmon directories
salmondirs <- list.files(salmondir, full.names = TRUE)
salmonfiles <- paste0(salmondirs, "/quant.sf")
names(salmonfiles) <- basename(salmondirs)

## Read transcript-to-gene mapping
tx2gene <- readRDS(tx2gene)

## Read Salmon abundances
txi <- tximport(files = salmonfiles, type = "salmon", txOut = TRUE)

## Filter out transcripts that are lowly expressed in all samples 
TPM <- as.data.frame(txi$abundance)
TPM$gene <- tx2gene$gene[match(rownames(TPM), tx2gene$tx)]
TPM$transcript <- tx2gene$tx[match(rownames(TPM), tx2gene$tx)]
tmp <- TPM %>% dplyr::group_by(gene) %>% dplyr::mutate_at(.funs = funs(perc = ./sum(.)), 
                                                          .vars = names(salmonfiles))
tmp[is.na(tmp)] <- 0
tmp <- as.data.frame(tmp)
tmp$gene <- NULL
tmp$maxfrac <- apply(tmp[, grep("perc", colnames(tmp))], 1, max)

keeptx <- tmp$transcript[tmp$maxfrac >= 0.05]
counts <- txi$counts[rownames(txi$counts) %in% keeptx, ]

## Read metadata and reorder in the same order as count matrix
metadata <- read.delim(metafile, header = TRUE, as.is = TRUE, sep = "\t")
metadata <- metadata[match(colnames(counts), metadata$ID), ]
conditions <- as.character(metadata$group)
names(conditions) <- colnames(counts)
print(conditions)

## Run DEXSeq
BPPARAM <- MulticoreParam(workers = 8)

## Generate DEXSeqDataSet
message("Generating DEXSeqDataSet...")
dxd0 <- DEXSeqDataSet(countData = round(counts), 
                      sampleData = data.frame(condition = conditions, 
                                              stringsAsFactors = FALSE),
                      design = ~ sample + exon + condition:exon, 
                      featureID = rownames(counts),
                      groupID = tx2gene$gene[match(rownames(counts), tx2gene$tx)])

## Define comparisons of interest (first value in vector, or the level 
## after "vs", will be the reference level)
comps <- list(Sprk1.vs.WT = c("WT", "Sprk1"))

dexseq_res <- lapply(comps, function(cm) {
  ## Subset DEXSeqDataSet
  message("Subsetting DEXSeqDataSet...")
  dxd <- dxd0[, dxd0$condition %in% cm]
  dxd$condition <- factor(dxd$condition, levels = cm)
  
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
  res <- res %>% dplyr::rename(gene_id = groupID, feature_id = featureID) %>%
    dplyr::left_join(tx2gene %>% dplyr::select(tx, symbol, gene_biotype) %>%
                       dplyr::distinct(), 
                     by = c("feature_id" = "tx")) %>%
    dplyr::select(gene_id, feature_id, symbol, gene_biotype, everything())
  
  res
})

## Save results
for (nm in names(dexseq_res)) {
  write.table(dexseq_res[[nm]] %>% dplyr::filter(padj <= 0.05) %>%
                dplyr::arrange(pvalue), 
              file = gsub("\\.rds", paste0("_", nm, "_signif.txt"), outrds), 
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}
saveRDS(list(dexseq_res = dexseq_res), file = outrds)

sessionInfo()
date()

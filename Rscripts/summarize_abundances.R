args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

## Summarize abundances on transcript- and gene-level

suppressPackageStartupMessages({
  library(tximport)
  library(dplyr)
})

print(topdir)
print(tx2gene)
print(outrds)

get_file_listing <- function(topdir, subdir, filename) {
  files <- list.files(paste0(topdir, "/", subdir), full.names = TRUE)
  nm <- basename(files)
  files <- paste0(files, "/", filename)
  names(files) <- nm
  files
}

import_salmon_files <- function(files, name, tx2gene) {
  txi <- tximport(files, type = "salmon", txOut = TRUE)
  txig <- summarizeToGene(txi, tx2gene = tx2gene[, c("tx", "gene")])
  txab <- 
    dplyr::full_join(
      txi$counts %>% as.data.frame() %>% 
        setNames(paste0(names(.), "__count__", name)) %>%
        tibble::rownames_to_column("tx"),
      txi$abundance %>% as.data.frame() %>%
        setNames(paste0(names(.), "__tpm__", name)) %>%
        tibble::rownames_to_column("tx"),
      by = "tx"
    )
  geneab <- 
    dplyr::full_join(
      txig$counts %>% as.data.frame() %>% 
        setNames(paste0(names(.), "__count__", name)) %>%
        tibble::rownames_to_column("gene"),
      txig$abundance %>% as.data.frame() %>%
        setNames(paste0(names(.), "__tpm__", name)) %>%
        tibble::rownames_to_column("gene"),
      by = "gene"
    )
  list(txab = txab, geneab = geneab)
}

## Read tx2gene
tx2gene <- readRDS(tx2gene)

## Initialize gene and transcript abundance lists
gene_abundances <- list()
tx_abundances <- list()

## Parse the subdirectories of topdir to see which abundances are available
(dirs <- list.files(topdir, full.names = FALSE))

## Salmon31
if ("salmon31" %in% dirs) {
  message("Salmon31")
  files <- get_file_listing(topdir, "salmon31", "quant.sf")
  print(files)
  abds <- import_salmon_files(files, "salmon31", tx2gene)
  tx_abundances[["salmon31"]] <- abds$txab
  gene_abundances[["salmon31"]] <- abds$geneab
}

## Salmonminimap2
if ("salmonminimap2" %in% dirs) {
  message("salmonminimap2")
  files <- get_file_listing(topdir, "salmonminimap2", "quant.sf")
  print(files)
  abds <- import_salmon_files(files, "salmonminimap2", tx2gene)
  tx_abundances[["salmonminimap2"]] <- abds$txab
  gene_abundances[["salmonminimap2"]] <- abds$geneab
}

## Salmonminimap2_p0.99
if ("salmonminimap2_p0.99" %in% dirs) {
  message("salmonminimap2_p0.99")
  files <- get_file_listing(topdir, "salmonminimap2_p0.99", "quant.sf")
  print(files)
  abds <- import_salmon_files(files, "salmonminimap2_p0.99", tx2gene)
  tx_abundances[["salmonminimap2_p0.99"]] <- abds$txab
  gene_abundances[["salmonminimap2_p0.99"]] <- abds$geneab
}

## wub 
if ("wubminimap2" %in% dirs) {
  message("wub")
  files <- get_file_listing(topdir, "wubminimap2", "bam_count_reads.tsv")
  print(files)
  tx_abundances[["wubminimap2"]] <- 
    Reduce(function(...) dplyr::full_join(..., by = "tx"),
           lapply(names(files), function(nm) {
             read.delim(files[nm], header = TRUE, as.is = TRUE) %>% 
               dplyr::select(Reference, Count) %>%
               setNames(c("tx", paste0(nm, "__count__wubminimap2")))
           })
    )
  gene_abundances[["wubminimap2"]] <- 
    dplyr::left_join(tx_abundances[["wubminimap2"]],
                     tx2gene[, c("tx", "gene")]) %>% 
    dplyr::select(-tx) %>%
    dplyr::group_by(gene) %>% dplyr::summarise_all(funs(sum)) %>%
    as.data.frame()
}

## featureCounts (primary) 
if ("featurecountsminimap2_primary" %in% dirs) {
  message("featureCounts (primary)")
  files <- get_file_listing(topdir, "featurecountsminimap2_primary", "featurecountsminimap2.txt")
  print(files)
  gene_abundances[["featurecountsminimap2primary"]] <- 
    Reduce(function(...) dplyr::full_join(..., by = "gene"),
           lapply(names(files), function(nm) {
             read.delim(files[nm], skip = 1, header = TRUE, as.is = TRUE) %>% 
               dplyr::select(-Chr, -Start, -End, -Strand, -Length) %>%
               setNames(c("gene", paste0(nm, "__count__featurecountsminimap2primary")))
           })
    )
}

## Remove version numbers from all genes and transcripts
tx_abundances <- lapply(tx_abundances, function(l) {
  l %>% dplyr::mutate(tx = gsub("\\.[0-9]+$", "", tx))
})
gene_abundances <- lapply(gene_abundances, function(l) {
  l %>% dplyr::mutate(gene = gsub("\\.[0-9]+$", "", gene))
})

## Merge all abundance estimates
tx_abundances <- Reduce(function(...) dplyr::full_join(..., by = "tx"),
                        tx_abundances)
gene_abundances <- Reduce(function(...) dplyr::full_join(..., by = "gene"),
                          gene_abundances)

## Set transcript/gene names as row names and replace NAs with zeros
rownames(tx_abundances) <- tx_abundances$tx
tx_abundances$tx <- NULL
rownames(gene_abundances) <- gene_abundances$gene
gene_abundances$gene <- NULL

tx_abundances[is.na(tx_abundances)] <- 0
gene_abundances[is.na(gene_abundances)] <- 0

saveRDS(list(tx_abundances = tx_abundances, 
             gene_abundances = gene_abundances),
        file = outrds)

date()
sessionInfo()


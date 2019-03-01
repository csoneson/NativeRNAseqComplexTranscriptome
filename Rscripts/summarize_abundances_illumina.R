args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

## Summarize abundances on transcript- and gene-level

suppressPackageStartupMessages({
  library(tximport)
  library(dplyr)
  library(rtracklayer)
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

## Salmon
if ("salmon31" %in% dirs) {
  message("Salmon")
  (files <- get_file_listing(topdir, "salmon31", "quant.sf"))
  abds <- import_salmon_files(files, "salmon", tx2gene)
  tx_abundances[["salmon"]] <- abds$txab
  gene_abundances[["salmon"]] <- abds$geneab
}

## StringTie
if ("stringtie" %in% dirs) {
  message("StringTie")
  files <- list.files(paste0(topdir, "/stringtie"), full.names = TRUE)
  nm <- basename(files)
  files <- paste0(files, "/", nm, "_stringtie.gtf")
  names(files) <- nm
  print(files)
  txi <- 
    Reduce(function(...) dplyr::full_join(..., by = c("tx", "gene")),
           lapply(names(files), function(nm) {
             message(nm)
             f <- rtracklayer::import(files[nm])
             f <- subset(f, type == "transcript")
             data.frame(tx = f$transcript_id, 
                        gene = f$gene_id,
                        tpm = as.numeric(as.character(f$TPM)),
                        stringsAsFactors = FALSE) %>%
               setNames(c("tx", "gene", paste0(nm, "__tpm__StringTie")))
           }))
  tx_abundances[["stringtie"]] <- txi %>% dplyr::select(-gene)
  gene_abundances[["stringtie"]] <- txi %>% dplyr::select(-tx) %>%
    dplyr::group_by(gene) %>% dplyr::summarize_all(funs(sum)) %>% as.data.frame()
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


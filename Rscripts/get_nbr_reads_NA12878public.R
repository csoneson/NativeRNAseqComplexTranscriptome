args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

## Catalog the number of nanopore reads

print(txomebamdir_p0.99)
print(outrds)

suppressPackageStartupMessages({
  library(ShortRead)
  library(GenomicAlignments)
  library(dplyr)
  library(ggplot2)
  library(data.table)
})

## ========================================================================== ##
## Help functions
## ========================================================================== ##
readBam <- function(bamfile) {
  bam <- readGAlignments(bamfile, use.names = TRUE,
                         param = ScanBamParam(tag = c("NM"),
                                              what = c("qname", "flag", "rname", 
                                                       "pos", "mapq")))
  mcols(bam)$tmp <- seq_len(length(bam)) %/% 5e5
  
  bamout <- do.call(c, lapply(unique(mcols(bam)$tmp), function(i) {
    tmp <- bam[mcols(bam)$tmp == i]
    ops <- GenomicAlignments::CIGAR_OPS
    wdths <- GenomicAlignments::explodeCigarOpLengths(cigar(tmp), ops = ops)
    keep.ops <- GenomicAlignments::explodeCigarOps(cigar(tmp), ops = ops)
    explodedcigars <- IRanges::CharacterList(relist(paste0(unlist(wdths), 
                                                           unlist(keep.ops)), wdths))
    for (opts in setdiff(GenomicAlignments::CIGAR_OPS, "=")) {
      mcols(tmp)[[paste0("nbr", opts)]] <- 
        sapply(explodedcigars, function(cg) sum(as.numeric(gsub(paste0(opts, "$"), "", cg)), na.rm = TRUE))
    }
    mcols(tmp)$readLength <- rowSums(as.matrix(mcols(tmp)[, c("nbrS", "nbrH", "nbrM", "nbrI")]))
    tmp
  }))
  
  bamout
}

makeReadDf <- function(bam) {
  tmp <- data.frame(bam %>% setNames(NULL), stringsAsFactors = FALSE) %>%
    dplyr::rename(read = qname,
                  nbrJunctions = njunc) %>%
    dplyr::select(-cigar) %>%
    dplyr::mutate(alignedLength = nbrM + nbrI) ## equivalent to readLength-nbrS-nbrH
  
  tmp2 <- as.data.frame(table(names(subset(bam, flag %in% c(0, 16)))))
  if (nrow(tmp2) == 0) tmp2 <- data.frame(Var1 = tmp$read[1], Freq = 0)
  tmp <- tmp %>% 
    dplyr::left_join(tmp2 %>% dplyr::rename(read = Var1, nbrPrimaryAlignments = Freq))
  
  tmp3 <- as.data.frame(table(names(subset(bam, flag %in% c(256, 272)))))
  if (nrow(tmp3) == 0) tmp3 <- data.frame(Var1 = tmp$read[1], Freq = 0)
  tmp <- tmp %>% 
    dplyr::left_join(tmp3 %>% dplyr::rename(read = Var1, nbrSecondaryAlignments = Freq))
  
  tmp %>% dplyr::mutate(nbrSecondaryAlignments = replace(nbrSecondaryAlignments, 
                                                         is.na(nbrSecondaryAlignments), 0))
}

makeSummaryList <- function(bam) {
  list(nAlignments = length(bam), 
       nAlignedReads = length(unique(names(bam))),
       nPrimaryAlignments = length(subset(bam, flag %in% c(0, 16))), 
       nReadsWithPrimaryAlignments = length(unique(names(subset(bam, flag %in% c(0, 16))))),
       nSecondaryAlignments = length(subset(bam, flag %in% c(256, 272))),
       nReadsWithSecondaryAlignments = length(unique(names(subset(bam, flag %in% c(256, 272))))))
}

## ========================================================================== ##
## BAM files, transcriptome alignment, -p 0.99
bamfilestx_p0.99 <- list.files(txomebamdir_p0.99, 
                               pattern = "_minimap_txome_p0.99.bam$", 
                               recursive = TRUE, full.names = TRUE)
names(bamfilestx_p0.99) <- gsub("_minimap_txome_p0.99.bam", "", 
                                basename(bamfilestx_p0.99))
bamfilestx_p0.99
txomebams_p0.99 <- lapply(bamfilestx_p0.99, function(f) {
  bam <- readBam(f)
  tmp <- makeReadDf(bam)
  
  ## Read transcript lengths
  bamheader <- scanBamHeader(f)[[1]]$targets
  txlengths <- data.frame(rname = names(bamheader),
                          txLength = bamheader,
                          stringsAsFactors = FALSE)
  
  c(makeSummaryList(bam),
    list(allAlignments = tmp %>% dplyr::left_join(txlengths, by = "rname")))
})
for (n in names(txomebams_p0.99)) {
  txomebams_p0.99[[n]][["allAlignments"]]$sample <- n
}

saveRDS(list(txomebams_p0.99 = txomebams_p0.99),
        file = outrds)
sessionInfo()
date()

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

## Catalog the number of nanopore reads

print(fastqdir)
print(genomebamdir)
print(txomebamdir)
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
                                              what = c("qname","flag", "rname", 
                                                       "pos", "mapq")))
  ops <- GenomicAlignments::CIGAR_OPS
  wdths <- GenomicAlignments::explodeCigarOpLengths(cigar(bam), ops = ops)
  keep.ops <- GenomicAlignments::explodeCigarOps(cigar(bam), ops = ops)
  explodedcigars <- IRanges::CharacterList(relist(paste0(unlist(wdths), 
                                                         unlist(keep.ops)), wdths))
  for (opts in setdiff(GenomicAlignments::CIGAR_OPS, "=")) {
    mcols(bam)[[paste0("nbr", opts)]] <- 
      sapply(explodedcigars, function(cg) sum(as.numeric(gsub(paste0(opts, "$"), "", cg)), na.rm = TRUE))
  }
  mcols(bam)$readLength <- rowSums(as.matrix(mcols(bam)[, c("nbrS", "nbrH", "nbrM", "nbrI")]))
  bam
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

  tmp4 <- as.data.frame(table(names(subset(bam, flag %in% c(2048, 2064)))))
  if (nrow(tmp4) == 0) tmp4 <- data.frame(Var1 = tmp$read[1], Freq = 0)
  tmp <- tmp %>% 
    dplyr::left_join(tmp4 %>% dplyr::rename(read = Var1, nbrSupplementaryAlignments = Freq))
  
  tmp %>% dplyr::mutate(nbrSecondaryAlignments = replace(nbrSecondaryAlignments, 
                                                         is.na(nbrSecondaryAlignments), 0),
                        nbrSupplementaryAlignments = replace(nbrSupplementaryAlignments, 
                                                             is.na(nbrSupplementaryAlignments), 0))
}

makeSummaryList <- function(bam) {
  list(nAlignments = length(bam), 
       nAlignedReads = length(unique(names(bam))),
       nPrimaryAlignments = length(subset(bam, flag %in% c(0, 16))), 
       nReadsWithPrimaryAlignments = length(unique(names(subset(bam, flag %in% c(0, 16))))),
       nSecondaryAlignments = length(subset(bam, flag %in% c(256, 272))),
       nReadsWithSecondaryAlignments = length(unique(names(subset(bam, flag %in% c(256, 272))))),
       nSupplementaryAlignments = length(subset(bam, flag %in% c(2048, 2064))),
       nReadsWithSupplementaryAlignments = length(unique(names(subset(bam, flag %in% c(2048, 2064)))))
  )
}

## ========================================================================== ##
## Read individual FASTQ files and get number of reads and read lengths
fastqfiles <- list.files(fastqdir, pattern = "(FASTQ|fastq)\\.gz$", full.names = TRUE)
names(fastqfiles) <- gsub("\\.(FASTQ|fastq).gz", "", basename(fastqfiles))
fastqfiles
fastqs <- lapply(fastqfiles, function(f) {
  fastq <- fread(paste0("zcat ", f), sep = "\n", header = FALSE)$V1
  seq_idxs <- seq(2, length(fastq), by = 4)
  list(nReads = length(seq_idxs), 
       reads = do.call(rbind, lapply(seq_idxs, function(i) {
         data.frame(read = gsub("^@", "", strsplit(fastq[i - 1], " ")[[1]][1]),
                    readLength = nchar(fastq[i]),
                    aveBaseQuality = mean(as.numeric(charToRaw(fastq[i + 2])) - 33),
                    stringsAsFactors = FALSE)
       })))
})

## BAM files, genome alignment
bamfiles <- list.files(genomebamdir, pattern = "_minimap_genome_s.bam$", 
                       recursive = TRUE, full.names = TRUE)
names(bamfiles) <- gsub("_minimap_genome_s.bam", "", basename(bamfiles))
bamfiles
genomebams <- lapply(bamfiles, function(f) {
  bam <- readBam(f)
  tmp <- makeReadDf(bam)
  c(makeSummaryList(bam), 
    list(allAlignments = tmp))
})
for (n in names(genomebams)) {
  genomebams[[n]][["allAlignments"]]$sample <- n
}

## BAM files, transcriptome alignment
bamfilestx <- list.files(txomebamdir, pattern = "_minimap_txome.bam$", 
                         recursive = TRUE, full.names = TRUE)
names(bamfilestx) <- gsub("_minimap_txome.bam", "", basename(bamfilestx))
bamfilestx
txomebams <- lapply(bamfilestx, function(f) {
  bam <- readBam(f)
  tmp <- makeReadDf(bam)
  
  ## Read transcript lengths
  bamheader <- scanBamHeader(f)[[1]]$targets
  txlengths <- data.frame(rname = names(bamheader),
                          txLength = bamheader,
                          stringsAsFactors = FALSE)
  
  c(makeSummaryList(bam), 
    list(allAlignments = tmp %>% dplyr::left_join(txlengths, by = "rname"),
         nTxCoveredByPrimary = length(unique(seqnames(subset(bam, flag %in% c(0, 16)))))))
})
for (n in names(txomebams)) {
  txomebams[[n]][["allAlignments"]]$sample <- n
}

## BAM files, transcriptome alignment, -p 0.99
bamfilestx_p0.99 <- list.files(txomebamdir_p0.99, pattern = "_minimap_txome_p0.99.bam$", 
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
    list(allAlignments = tmp %>% dplyr::left_join(txlengths, by = "rname"),
         nTxCoveredByPrimary = length(unique(seqnames(subset(bam, flag %in% c(0, 16)))))))
})
for (n in names(txomebams_p0.99)) {
  txomebams_p0.99[[n]][["allAlignments"]]$sample <- n
}

nReadTable <- 
  data.frame(sample = names(fastqs), 
             nReads = sapply(fastqs, function(w) w$nReads),
             stringsAsFactors = FALSE) %>%
  dplyr::full_join(data.frame(
    sample = names(genomebams), 
    nAlignedReadsGenome = sapply(genomebams, function(w) w$nAlignedReads),
    nAlignmentsGenome = sapply(genomebams, function(w) w$nAlignments),
    nPrimaryAlignmentsGenome = sapply(genomebams, function(w) w$nPrimaryAlignments),
    nReadsWithPrimaryAlignmentsGenome = sapply(genomebams, function(w) w$nReadsWithPrimaryAlignments),
    nSecondaryAlignmentsGenome = sapply(genomebams, function(w) w$nSecondaryAlignments),
    nReadsWithSecondaryAlignmentsGenome = sapply(genomebams, function(w) w$nReadsWithSecondaryAlignments),
    nSupplementaryAlignmentsGenome = sapply(genomebams, function(w) w$nSupplementaryAlignments),
    nReadsWithSupplementaryAlignmentsGenome = sapply(genomebams, function(w) w$nReadsWithSupplementaryAlignments),
    stringsAsFactors = FALSE)) %>%
  dplyr::full_join(data.frame(
    sample = names(txomebams), 
    nAlignedReadsTxome = sapply(txomebams, function(w) w$nAlignedReads),
    nAlignmentsTxome = sapply(txomebams, function(w) w$nAlignments),
    nTxCoveredByPrimary = sapply(txomebams, function(w) w$nTxCoveredByPrimary),
    nPrimaryAlignmentsTxome = sapply(txomebams, function(w) w$nPrimaryAlignments),
    nReadsWithPrimaryAlignmentsTxome = sapply(txomebams, function(w) w$nReadsWithPrimaryAlignments),
    nSecondaryAlignmentsTxome = sapply(txomebams, function(w) w$nSecondaryAlignments),
    nReadsWithSecondaryAlignmentsTxome = sapply(txomebams, function(w) w$nReadsWithSecondaryAlignments),
    nSupplementaryAlignmentsTxome = sapply(txomebams, function(w) w$nSupplementaryAlignments),
    nReadsWithSupplementaryAlignmentsTxome = sapply(txomebams, function(w) w$nReadsWithSupplementaryAlignments),
    stringsAsFactors = FALSE)) %>%
  dplyr::full_join(data.frame(
    sample = names(txomebams_p0.99), 
    nAlignedReadsTxome_p0.99 = sapply(txomebams_p0.99, function(w) w$nAlignedReads),
    nAlignmentsTxome_p0.99 = sapply(txomebams_p0.99, function(w) w$nAlignments),
    nTxCoveredByPrimary_p0.99 = sapply(txomebams_p0.99, function(w) w$nTxCoveredByPrimary),
    nPrimaryAlignmentsTxome_p0.99 = sapply(txomebams_p0.99, function(w) w$nPrimaryAlignments),
    nReadsWithPrimaryAlignmentsTxome_p0.99 = sapply(txomebams_p0.99, function(w) w$nReadsWithPrimaryAlignments),
    nSecondaryAlignmentsTxome_p0.99 = sapply(txomebams_p0.99, function(w) w$nSecondaryAlignments),
    nReadsWithSecondaryAlignmentsTxome_p0.99 = sapply(txomebams_p0.99, function(w) w$nReadsWithSecondaryAlignments),
    nSupplementaryAlignmentsTxome_p0.99 = sapply(txomebams_p0.99, function(w) w$nSupplementaryAlignments),
    nReadsWithSupplementaryAlignmentsTxome_p0.99 = sapply(txomebams_p0.99, function(w) w$nReadsWithSupplementaryAlignments),
    stringsAsFactors = FALSE))

write.table(nReadTable, file = gsub("rds$", "txt", outrds), 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

saveRDS(list(fastqs = fastqs, genomebams = genomebams, txomebams = txomebams,
             txomebams_p0.99 = txomebams_p0.99, nReadTable = nReadTable),
        file = outrds)
sessionInfo()
date()

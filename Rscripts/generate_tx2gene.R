args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

## Parse a transcript fasta file and optionally a gtf file from Ensembl and
## extract information about the transcripts.

suppressPackageStartupMessages({
  library(Biostrings)
  library(rtracklayer)
})

print(transcriptome)
print(gtf)
print(outrds)

## Read input files
transcriptome <- readDNAStringSet(transcriptome)
if (!is.null(gtf)) gtf0 <- import(gtf)

## Extract information from fasta file identifiers
tx2gene <- data.frame(t(sapply(as.character(names(transcriptome)), function(nm) {
  a <- strsplit(nm, " ")[[1]]
  tx <- a[1]
  gene <- gsub("^gene:", "", a[grep("^gene:", a)])
  symbol <- gsub("^gene_symbol:", "", a[grep("^gene_symbol:", a)])
  gene_biotype <- gsub("^gene_biotype:", "", a[grep("^gene_biotype:", a)])
  tx_biotype <- gsub("^transcript_biotype:", "", a[grep("^transcript_biotype:", a)])
  position <- gsub("chromosome:|scaffold:", "", a[grep("^chromosome:|^scaffold:", a)])
  c(tx = ifelse(length(tx) != 0, tx, NA), 
    gene = ifelse(length(gene) != 0, gene, NA),
    symbol = ifelse(length(symbol) != 0, symbol, NA),
    gene_biotype = ifelse(length(gene_biotype) != 0, gene_biotype, NA),
    tx_biotype = ifelse(length(tx_biotype) != 0, tx_biotype, NA),
    chromosome = ifelse(length(position) != 0, strsplit(position, ":")[[1]][2], NA),
    start = ifelse(length(position) != 0, strsplit(position, ":")[[1]][3], NA),
    end = ifelse(length(position) != 0, strsplit(position, ":")[[1]][4], NA),
    strand = ifelse(length(position) != 0, strsplit(position, ":")[[1]][5], NA))
})), stringsAsFactors = FALSE)
rownames(tx2gene) <- NULL

## Add information about transcript length
tx2gene$txlength <- width(transcriptome)

if (!is.null(gtf)) {
  ## Add information from gtf file
  gtf <- subset(gtf0, type == "transcript")
  gtf$transcript_id_with_version <- paste0(gtf$transcript_id, ".", gtf$transcript_version)
  ## Keep only transcripts that are not present in the fasta file
  gtf <- subset(gtf, !(transcript_id_with_version %in% tx2gene$tx))
  gtfout <- data.frame(tx = gtf$transcript_id_with_version, 
                       gene = paste0(gtf$gene_id, ".", gtf$gene_version),
                       symbol = gtf$gene_name,
                       gene_biotype = gtf$gene_biotype,
                       tx_biotype = gtf$transcript_biotype,
                       chromosome = seqnames(gtf),
                       start = start(gtf),
                       end = end(gtf),
                       strand = strand(gtf),
                       stringsAsFactors = FALSE)
  
  ## Get transcript length by adding up exon lengths for each transcript
  gtfe <- subset(gtf0, type == "exon" & 
                   paste0(transcript_id, ".", transcript_version) %in% gtfout$tx)
  gtfout$txlength <- sapply(gtfout$tx, function(m) {
    sum(width(subset(gtfe, paste0(transcript_id, ".", transcript_version) == m)))
  })
  
  ## Merge information from fasta and gtf
  stopifnot(all(colnames(tx2gene) == colnames(gtfout)))
  tx2gene <- rbind(tx2gene, gtfout)
}

saveRDS(tx2gene, file = outrds)

sessionInfo()
date()

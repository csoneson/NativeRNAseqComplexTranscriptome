args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(bamfile)
print(fastqfile)
print(outrds)

suppressPackageStartupMessages({
  library(GenomicAlignments)
  library(data.table)
  library(ggplot2)
})

sbp <- ScanBamParam(tag = "NM")

## Read bam file and corresponding fastq file
bam <- readGAlignments(bamfile, use.names = TRUE, param = sbp)
fastq <- fread(paste0("gunzip -c ", fastqfile), sep = "\n", header = FALSE)

## Split rows in fastq file
id_idxs <- seq(1, nrow(fastq), by = 4)   ## read ids
q_idxs <- seq(4, nrow(fastq), by = 4)   ## quality scores
seq_idxs <- seq(2, nrow(fastq), by = 4)   ## sequence

## Calculate average quality per read
mean_q <- sapply(fastq$V1[q_idxs], function(u) mean(as.numeric(charToRaw(u)) - 33))
names(mean_q) <- gsub("@", "", sapply(strsplit(fastq$V1[id_idxs]," "), .subset, 1))

## Get read lengths
seqlen <- nchar(fastq$V1[seq_idxs])

## Number of soft clipped bases
clip_length <- sapply(explodeCigarOpLengths(cigar(bam), ops = "S"), sum)

pdf(gsub("rds$", "pdf", outrds), w = 10, h = 10)
## mapped length vs number of mismatches
smoothScatter(log10(qwidth(bam) - clip_length), log10(mcols(bam)$NM), 
               xlab = "mapped length (log10)", ylab = "number of mismatches (log10)")

## percent mismatches
pct_mm <- mcols(bam)$NM/(qwidth(bam) - clip_length)
hist(pct_mm[pct_mm >= 0 & pct_mm <= 1], 500, main = "percent mismatches")

## read length
hist(log10(seqlen), 500, main = "read length (log10)")

## mapped length
hist(log10(qwidth(bam) - clip_length), 500, main = "mapped length (log10)")

## read length vs average quality score
smoothScatter(log10(seqlen), mean_q, xlab = "read length (log10)", 
              ylab = "average quality score")

df <- data.frame(sequence_length = seqlen, mean_quality = mean_q)

p <- ggplot(df) + geom_hex(aes(sequence_length, mean_quality), bins = 100) +
         scale_x_log10() +
         scale_fill_gradientn("", colours = rev(rainbow(10, end = 4/6))) +
         theme(text = element_text(size=20))
print(p)

## average quality score
hist(mean_q, 100, main = "average quality score")

## average quality score vs percent mismatches
m <- match(names(bam), names(mean_q))
smoothScatter(mean_q[m], pct_mm, pch = ".", ylim = c(0, 0.4), 
              xlab = "average quality score", 
              ylab = "percent mismatches")

df <- data.frame(mean_quality = mean_q[m], percent_mismatches=pct_mm*100)

p <- ggplot(df) + geom_hex(aes(mean_quality, percent_mismatches), bins = 100) +
         scale_x_log10() +
         scale_fill_gradientn("", colours = rev(rainbow(10, end = 4/6))) +
         theme(text = element_text(size=20))
print(p)

dev.off()

saveRDS(NULL, outrds)
sessionInfo()
date()

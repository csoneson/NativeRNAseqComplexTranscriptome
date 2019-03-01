args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(refgtf)
print(bam)
print(paired)
print(outrds)

suppressPackageStartupMessages({
  library(GenomicAlignments)
  library(SGSeq)
  library(GenomicFeatures)
})

## Extract annotated junctions
txdb <- GenomicFeatures::makeTxDbFromGFF(refgtf)
txf <- SGSeq::convertToTxFeatures(txdb)
txf <- txf[SGSeq::type(txf) == "J"]
start(txf) <- start(txf) + 1
end(txf) <- end(txf) - 1

## Get junctions from alignments
if (!paired) {
  x <- GenomicAlignments::readGAlignments(bam, use.names = TRUE)
} else {
  x <- GenomicAlignments::readGAlignmentPairs(bam, use.names = TRUE)
}
juncSum <- GenomicAlignments::summarizeJunctions(x)

## To compare with later: number of exact matches of observed junctions among
## the annotated ones.
table(ranges(juncSum) %in% ranges(txf))

## Find closest annotated junction to each observed junction
txf_seqnames_char <- as.character(seqnames(txf))
juncsum_seqnames_char <- as.character(seqnames(juncSum))
txf_starts <- start(txf)
txf_ends <- end(txf)
juncsum_starts <- start(juncSum)
juncsum_ends <- end(juncSum)
juncdists <- do.call(rbind, lapply(seq_along(juncSum), function(i) {
  if (i%%1000 == 0) print(i)
  idx <- which(txf_seqnames_char == juncsum_seqnames_char[i])
  dists <- abs(juncsum_starts[i] - txf_starts[idx]) + 
    abs(juncsum_ends[i] - txf_ends[idx])
  tryCatch({
    data.frame(juncSumIdx = i, whichMin = idx[which.min(dists)], minDist = min(dists))
  }, error = function(e) {
    data.frame(juncSumIdx = i, whichMin = NA_integer_, minDist = Inf)
  })
}))

rtmp <- ranges(txf)
l <- length(rtmp)
rtmp <- c(rtmp, IRanges(start = 0, end = 0))
jtmp <- juncdists
jtmp$whichMin[is.na(jtmp$whichMin)] <- l + 1
mcols(juncSum)$closestRefJunc <- rtmp[jtmp$whichMin]
mcols(juncSum)$distToClosestRefJunc <- juncdists$minDist

saveRDS(juncSum, file = outrds)
date()
sessionInfo()

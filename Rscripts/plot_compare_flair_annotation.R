args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(annotatedGtf)
print(flairGtf)
print(dataSet)
print(outDir)
print(seed)
print(nPlots)
print(flairTmap)  ## gffcompare classification
print(sqanticlass)  ## SQANTI classification
print(outrds)

suppressPackageStartupMessages({
  library(S4Vectors)
  library(rtracklayer)
  library(Gviz)
  library(dplyr)
})

## Read gtf files
annotated <- rtracklayer::import(annotatedGtf)
idx <- match(c("transcript_id", "gene_id", "exon_id"), 
             colnames(mcols(annotated)))
colnames(mcols(annotated))[idx] <- c("transcript", "gene", "exon")
S4Vectors::mcols(annotated)$symbol <- mcols(annotated)$transcript
annotated <- subset(annotated, type == "exon")

flair <- rtracklayer::import(flairGtf)
colnames(mcols(flair))[match("transcript_id", colnames(mcols(flair)))] <- 
  "transcript"
S4Vectors::mcols(flair)$symbol <- mcols(flair)$transcript
flair <- subset(flair, type == "exon")

## Read tmap file to sample transcript pairs from and sample nPlots pairs
## Also read SQANTI classification info to add on top
tmap <- read.delim(flairTmap, header = TRUE, as.is = TRUE)
if (sqanticlass != "") {
  sqanti <- read.delim(sqanticlass, header = TRUE, as.is = TRUE)
}
if (sqanticlass != "") {
  tmap <- tmap %>% dplyr::left_join(
    sqanti %>% dplyr::select(isoform, associated_transcript, structural_category),
    by = c("qry_id" = "isoform")
  ) %>%
    dplyr::filter(ref_id == associated_transcript | 
                    !(class_code %in% c("c", "=")))
} else {
  tmap$structural_category <- ""
}
set.seed(seed)
tmap <- tmap[sample(seq_len(nrow(tmap)), size = nPlots, replace = FALSE), ]

## Define help function
compareFlairAnnotatedTx <- function(dataSet, annotatedGtf, flairGtf, 
                                    annotatedTx, flairTx, classCode, 
                                    structCat, outDir) {
  options(ucscChromosomeNames = FALSE)
  
  ## Subset gtf files to selected transcripts
  atmp <- subset(annotatedGtf, transcript == annotatedTx)
  ftmp <- subset(flairGtf, transcript == flairTx)
  
  ## Get width of each transcript
  awd <- sum(width(atmp))
  fwd <- sum(width(ftmp))
  
  show_chr <- unique(c(seqnames(atmp), seqnames(ftmp)))
  stopifnot(length(show_chr) == 1)
  
  min_coord <- min(
    min(start(atmp)) - 0.2*(max(end(atmp)) - min(start(atmp))),
    min(start(ftmp)) - 0.2*(max(end(ftmp)) - min(start(ftmp)))
  )
  max_coord <- max(
    max(end(atmp)) + 0.05*(max(end(atmp)) - min(start(atmp))),
    max(end(ftmp)) + 0.05*(max(end(ftmp)) - min(start(ftmp)))
  )
  
  agrtr <- Gviz::GeneRegionTrack(atmp, showId = FALSE, col = NULL,
                                 fill = "blue", name = "annotated",
                                 background.title = "transparent",
                                 col.title = "black", min.height = 15)
  fgrtr <- Gviz::GeneRegionTrack(ftmp, showId = FALSE, col = NULL,
                                 fill = "orange", name = "flair", 
                                 background.title = "transparent",
                                 col.title = "black", min.height = 15)
  
  gtr <- Gviz::GenomeAxisTrack()
  
  tracks <- c(gtr, agrtr, fgrtr)
  
  pdf(paste0(outDir, "/", dataSet, "__flair_comparison__", 
             gsub(":", "_", annotatedTx), "__", 
             gsub(":", "_", flairTx), "__", classCode, "__", 
             structCat, ".pdf"),
      height = 3, width = 10)
  Gviz::plotTracks(tracks, chromosome = show_chr,
                   from = min_coord, to = max_coord, 
                   main = paste0("annotated: ", annotatedTx, " (", awd,
                                 " nt); flair: ", flairTx, " (", fwd, 
                                 " nt); class code: ", classCode, 
                                 "; structural category: ", structCat),
                   cex.main = 0.7,
                   min.width = 0, min.distance = 0, collapse = FALSE)
  dev.off()
}

options(ucscChromosomeNames = FALSE)
for (i in seq_len(nrow(tmap))) {
  compareFlairAnnotatedTx(dataSet = dataSet, annotatedGtf = annotated, 
                          flairGtf = flair, 
                          annotatedTx = tmap$ref_id[i], 
                          flairTx = tmap$qry_id[i], 
                          classCode = tmap$class_code[i],
                          structCat = tmap$structural_category[i], 
                          outDir = outDir)
}

saveRDS(tmap, file = outrds)
date()
sessionInfo()

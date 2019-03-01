args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

## Find the distance between supplementary alignments and the corresponding
## primary alignment

print(genomebamdir)
print(outrds)

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomicAlignments)
  library(dplyr)
  library(ggplot2)
})

bamfiles <- list.files(genomebamdir, pattern = "_minimap_genome_s.bam$", 
                       recursive = TRUE, full.names = TRUE)
names(bamfiles) <- gsub("_minimap_genome_s.bam", "", basename(bamfiles))
bamfiles

primary_supplementary_distances <- lapply(bamfiles, function(bam) {
  bam <- readGAlignments(bam, use.names = TRUE, 
                         param = ScanBamParam(tag=c("NM"),
                                              what=c("qname","flag", "rname", 
                                                     "pos", "mapq")))
  bam <- subset(bam, flag %in% c(0, 16, 2048, 2064))
  
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
  
  supplementary_alignments <- subset(bam, flag %in% c(2048, 2064))
  primary_alignments <- subset(bam, flag %in% c(0, 16))
  primary_alignments <- primary_alignments[names(primary_alignments) %in% 
                                             names(supplementary_alignments)]
  
  do.call(rbind, lapply(seq_along(supplementary_alignments), function(i) {
    nm <- names(supplementary_alignments)[i]
    primary <- as(primary_alignments[nm], "GRanges")
    supplementary <- as(supplementary_alignments[i], "GRanges")
    if (all(seqnames(primary) == seqnames(supplementary))) {
      ## Same chromosome
      if (all(strand(primary) == strand(supplementary))) {
        strands <- "same chromosome, same strand"
      } else {
        strands <- "same chromosome, different strand"
      }
      distn <- GenomicRanges::distance(primary, supplementary, ignore.strand = TRUE)
      ovlap <- sum(width(GenomicRanges::intersect(as(primary_alignments[nm], "GRangesList"),
                                                  as(supplementary_alignments[i], "GRangesList"),
                                                  ignore.strand = TRUE)))
    } else {
      strands <- "different chromosomes"
      distn <- NA
      ovlap <- NA
    }
    nbrMPrimary <- mcols(primary)[["nbrM"]]
    nbrMSupplementary <- mcols(supplementary)[["nbrM"]]
    nbrMSuppVsPrim <- nbrMSupplementary/nbrMPrimary
    data.frame(read = nm, strands = strands, distn = distn, ovlap = ovlap, 
               nbrMPrimary = nbrMPrimary, nbrMSupplementary = nbrMSupplementary,
               nbrMSuppVsPrim = nbrMSuppVsPrim, stringsAsFactors = FALSE)
  }))
})

for (nm in names(primary_supplementary_distances)) {
  primary_supplementary_distances[[nm]]$sample <- nm
}

primary_supplementary_distances <- do.call(rbind, primary_supplementary_distances)

pdf(gsub("rds$", "pdf", outrds), width = 8, height = 7)
ggplot(primary_supplementary_distances %>% 
         dplyr::mutate(dist0 = "") %>%
         dplyr::mutate(dist0 = replace(dist0, distn == 0, ", distance = 0"),
                       dist0 = replace(dist0, distn > 0, ", distance > 0")) %>% 
         dplyr::group_by(sample, strands, dist0) %>% dplyr::tally() %>%
         dplyr::mutate(strands_dist0 = paste0(strands, dist0)),
       aes(x = sample, y = n, fill = strands_dist0)) + 
  theme_bw() + geom_bar(stat = "identity") + 
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  ylab("Number of primary/supplementary alignment pairs") + xlab("") + 
  guides(fill = guide_legend(nrow = 5, byrow = TRUE)) + 
  scale_fill_manual(values = c("#DC050C", "#1965B0", "#4EB265", "#777777", "#B17BA6"),
                    name = "")
dev.off()

saveRDS(primary_supplementary_distances, file = outrds)
date()
sessionInfo()



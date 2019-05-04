args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(GenomicAlignments)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(stringr)
})
source("manuscript_results_figures/Rscripts/remap_sample_names.R")

conditions <- strsplit(conditions, ",")[[1]]

print(conditions)
print(genomebamdir)
print(outrds)

readBam <- function(bamfile) {
  bf <- BamFile(bamfile, yieldSize = 5e6)
  bam <- readGAlignments(bf, use.names = TRUE,
                         param = ScanBamParam(tag = c("NM"),
                                              what = c("qname", "flag", "rname", 
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
    dplyr::filter(flag %in% c(83, 99, 339, 355)) %>% ## mapped in proper pair, first in pair
    dplyr::mutate(alignedLength = nbrM + nbrI) ## equivalent to readLength-nbrS-nbrH
  
  tmp2 <- as.data.frame(table(names(subset(bam, flag %in% c(83, 99)))))
  if (nrow(tmp2) == 0) tmp2 <- data.frame(Var1 = tmp$read[1], Freq = 0)
  tmp <- tmp %>% 
    dplyr::left_join(tmp2 %>% dplyr::rename(read = Var1, nbrPrimaryAlignments = Freq))
  
  tmp3 <- as.data.frame(table(names(subset(bam, flag %in% c(339, 355)))))
  if (nrow(tmp3) == 0) tmp3 <- data.frame(Var1 = tmp$read[1], Freq = 0)
  tmp <- tmp %>% 
    dplyr::left_join(tmp3 %>% dplyr::rename(read = Var1, nbrSecondaryAlignments = Freq))
  
  tmp %>% dplyr::mutate(nbrSecondaryAlignments = replace(nbrSecondaryAlignments, 
                                                         is.na(nbrSecondaryAlignments), 0))
}

bamfiles <- list.files(genomebamdir, pattern = "_Aligned.sortedByCoord.out.bam$", 
                       recursive = TRUE, full.names = TRUE)
names(bamfiles) <- gsub("_Aligned.sortedByCoord.out.bam", "", basename(bamfiles))
bamfiles <- bamfiles[names(bamfiles) %in% 
                       sample_annotation$sample_orig[sample_annotation$condition %in% conditions]]
bamfiles
genomebams <- lapply(bamfiles, function(f) {
  bam <- readBam(f)
  tmp <- makeReadDf(bam)
  list(allAlignments = tmp)
})
for (n in names(genomebams)) {
  genomebams[[n]][["allAlignments"]]$sample <- n
}

genomebams <- do.call(dplyr::bind_rows, lapply(names(genomebams), function(nm) {
  genomebams[[nm]]$allAlignments %>% 
    dplyr::mutate(sample = remap[nm], dataset = "Illumina", 
                  rtype = "Reads aligning to the genome") %>%
    dplyr::mutate(fracM = nbrM/readLength,
                  fracS = nbrS/readLength,
                  fracI = nbrI/readLength)
}))

djunc <- genomebams %>% dplyr::filter(flag %in% c(83, 99)) %>%
  dplyr::mutate(nbrJunctions = replace(nbrJunctions, nbrJunctions > 20, ">20")) %>%
  dplyr::mutate(nbrJunctions = factor(nbrJunctions, 
                                      levels = c(as.character(0:20), ">20"))) %>%
  dplyr::group_by(sample, dataset, nbrJunctions) %>%
  dplyr::tally() %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(n = n/sum(n))

png(gsub("\\.rds$", "_njunc_distribution_bysample.png", outrds),
    width = 12, height = 8.5, unit = "in", res = 400)
print(ggplot(djunc,
             aes(x = nbrJunctions, y = n, fill = dataset, color = dataset)) +
        geom_bar(stat = "identity", position = "dodge") +
        theme_bw() + facet_wrap(~ sample) + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
              legend.position = "none") +
        scale_color_manual(values = ds_colors, name = "") + 
        scale_fill_manual(values = ds_colors, name = "") + 
        scale_y_continuous(expand = c(0, 0, 0.05, 0)) + 
        xlab("Number of spanned junctions") + ylab("Fraction of reads"))
dev.off()

png(gsub("\\.rds$", "_njunc_distribution.png", outrds),
    width = 8, height = 6, unit = "in", res = 400)
print(ggplot(djunc, aes(x = nbrJunctions, y = n), fill = "#B3B3B3", color = "#B3B3B3") +
        geom_bar(stat = "identity", position = "dodge") +
        theme_bw() + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
              legend.position = "none") +
        scale_y_continuous(expand = c(0, 0, 0.05, 0)) + 
        xlab("Number of spanned junctions") + ylab("Fraction of reads"))
dev.off()

saveRDS(djunc, file = outrds)
date()
sessionInfo()

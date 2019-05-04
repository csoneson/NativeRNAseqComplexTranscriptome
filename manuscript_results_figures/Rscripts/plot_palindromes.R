args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(stringr)
  library(Biostrings)
  library(data.table)
  library(parallel)
  library(BiocParallel)
})
source("manuscript_results_figures/Rscripts/remap_sample_names.R")

datasets <- strsplit(datasets, ",")[[1]]
names(datasets) <- datasets
conditions <- strsplit(conditions, ",")[[1]]

print(datasets)
print(conditions)
print(nreads)
print(ncores)
print(outrds)

## Information about the reads
readInfo <- lapply(datasets, function(ds) {
  rd <- readRDS(paste0(ds, "/output/", ds, "_nbr_reads.rds"))
  
  rdf <- rd$fastqs
  rdf <- rdf[!grepl("_orig", names(rdf))]
  rdf <- rdf[names(rdf) %in% 
               sample_annotation$sample_orig[sample_annotation$condition %in% conditions]]
  rdf <- do.call(dplyr::bind_rows, lapply(names(rdf), function(nm) {
    rdf[[nm]]$reads %>% 
      dplyr::mutate(sample = remap[nm], dataset = remapds[ds], rtype = "All reads")
  }))
  
  rdg <- rd$genomebams
  rdg <- rdg[!grepl("_orig", names(rdg))]
  rdg <- rdg[names(rdg) %in% 
               sample_annotation$sample_orig[sample_annotation$condition %in% conditions]]
  rdg <- do.call(dplyr::bind_rows, lapply(names(rdg), function(nm) {
    rdg[[nm]]$allAlignments %>% 
      dplyr::mutate(sample = remap[nm], dataset = remapds[ds], 
                    rtype = "Reads aligning to the genome") %>%
      dplyr::mutate(fracM = nbrM/readLength,
                    fracS = nbrS/readLength,
                    fracI = nbrI/readLength,
                    accuracy = (nbrM - NM + nbrI + nbrD)/(nbrM + nbrI + nbrD))
  }))
  
  list(rdf = rdf, rdg = rdg)
})

## All reads
allReads <- do.call(dplyr::bind_rows, lapply(readInfo, function(x) x$rdf))

## Genome alignments
gAlign <- do.call(dplyr::bind_rows, lapply(readInfo, function(x) x$rdg))

## Add information about genome alignments to the data frame with all reads
allReads <- allReads %>%
  dplyr::left_join(gAlign %>% dplyr::filter(flag %in% c(0, 16)) %>%
                     dplyr::select(read, sample, dataset, NM, nbrM, nbrS, 
                                   nbrD, nbrI, alignedLength, 
                                   nbrSecondaryAlignments, 
                                   nbrSupplementaryAlignments, fracM, fracS),
                   by = c("sample", "dataset", "read")) %>%
  dplyr::mutate(aligned = ifelse(is.na(alignedLength), "Not aligned", 
                                 "Aligned"))

## Get palindrome information for a subset of reads for each sample
fastqs <- do.call(dplyr::bind_rows, lapply(datasets, function(ds) {
  fastqdir <- paste0(ds, "/FASTQ", 
                     ifelse(ds %in% c("RNA001", "HEK293RNA"), "dna", ""))
  fastqfiles <- list.files(fastqdir, pattern = "(FASTQ|fastq)\\.gz$", 
                           full.names = TRUE)
  names(fastqfiles) <- gsub("\\.(FASTQ|fastq).gz", "", basename(fastqfiles))
  fastqfiles <- fastqfiles[!grepl("_orig", names(fastqfiles))]
  fastqfiles <- fastqfiles[names(fastqfiles) %in% 
                             sample_annotation$sample_orig[sample_annotation$condition %in% conditions]]
  fastqfiles
  
  do.call(dplyr::bind_rows, lapply(names(fastqfiles), function(nm) {
    f <- fastqfiles[nm]
    fastq <- fread(paste0("zcat ", f), sep = "\n", header = FALSE)$V1
    seq_idxs <- seq(2, length(fastq), by = 4)
    reads <- Biostrings::DNAStringSet(x = fastq[seq_idxs])
    names(reads) <- sapply(strsplit(fastq[seq_idxs - 1], " "), .subset, 1)
    print(nm)
    print(table(width(reads) < 30000))
    reads <- reads[width(reads) < 30000]
    set.seed(1)
    reads <- reads[sample(seq_along(reads), nreads, replace = FALSE)]
    readlengths <- width(reads)
    readnames <- names(reads)
    do.call(dplyr::bind_rows, bplapply(seq_along(reads), function(i) {
      fp <- Biostrings::findPalindromes(reads[[i]], max.looplength = readlengths[i], 
                                        min.armlength = 10)
      armlength <- Biostrings::palindromeArmLength(fp)
      totlength <- width(fp)
      w <- which.max(armlength)
      
      ## Get the total length of the read that is covered by palindromic sequence
      dnatmp <- reduce(IRanges(start = start(fp), 
                               width = armlength))
      totalPalindromeLength <- 2 * sum(width(dnatmp))
      tryCatch({
        data.frame(
          read = gsub("^@", "", readnames[i]),
          readLength = readlengths[i],
          maxArmLengthPalindrome = armlength[w],
          totalPalindromeLength = totlength[w],
          longestPalindromeArm = Biostrings::palindromeLeftArm(fp[w]),
          stringsAsFactors = FALSE
        )}, error = function(e) {
          data.frame(read = gsub("^@", "", readnames[i]),
                     readLength = readlengths[i],
                     maxArmLengthPalindrome = NA,
                     totalPalindromeLength = NA, 
                     longestPalindromeArm = "",
                     stringsAsFactors = FALSE)
        })
    }, BPPARAM = MulticoreParam(workers = ncores))) %>% 
      dplyr::mutate(sample = remap[nm], 
                    dataset = remapds[ds])
  }))
}))

df <- fastqs %>%
  dplyr::left_join(allReads, 
                   by = c("dataset", "sample", "read", "readLength")) %>%
  dplyr::mutate(
    maxArmLengthPalindrome = replace(maxArmLengthPalindrome, 
                                     is.na(maxArmLengthPalindrome), 0),
    totalPalindromeLength = replace(totalPalindromeLength, 
                                    is.na(totalPalindromeLength), 0))

## Read information about primary/supplementary alignments
dfprimsupp <- do.call(dplyr::bind_rows, lapply(datasets, function(ds) {
  readRDS(paste0(ds, "/output/", ds, 
                 "_primary_supplementary_alignments_distances.rds")) %>% 
    dplyr::mutate(sample = remap[sample]) %>% 
    dplyr::mutate(dataset = remapds[ds])
})) %>%
  dplyr::left_join(sample_annotation %>% dplyr::select(sample_remap, condition),
                   by = c("sample" = "sample_remap")) %>%
  dplyr::filter(condition %in% conditions) %>%
  dplyr::select(-condition) %>% 
  dplyr::mutate(dist0 = "") %>%
  dplyr::mutate(dist0 = replace(dist0, distn == 0, ", overlapping"),
                dist0 = replace(dist0, distn > 0, ", non-overlapping")) %>%
  dplyr::mutate(strands_dist0 = paste0(strands, dist0))

df <- df %>%
  dplyr::left_join(dfprimsupp, 
                   by = c("dataset", "sample", "read")) %>%
  dplyr::mutate(strands_dist0 = replace(strands_dist0, is.na(strands_dist0),
                                        "no supplementary alignment")) %>%
  dplyr::mutate(relMaxArmLengthPalindrome = maxArmLengthPalindrome/readLength)

catcols <- c("no supplementary alignment" = "lightyellow", 
             "different chromosomes" = "#E8601C", 
             "same chromosome, different strand, non-overlapping" = "#7BAFDE",
             "same chromosome, different strand, overlapping" = "#90C987", 
             "same chromosome, same strand, non-overlapping" = "#777777", 
             "same chromosome, same strand, overlapping" = "#B17BA6")

df <- df %>%
  dplyr::mutate(strands_dist0 = factor(strands_dist0, 
                                       levels = names(catcols)))

png(gsub("\\.rds$", "_aligned_vs_relpalindromelength.png", outrds),
    width = 10, height = 5, unit = "in", res = 400)
ggplot(df, aes(x = aligned, y = relMaxArmLengthPalindrome)) + 
  geom_violin() + geom_boxplot(width = 0.1) + 
  theme_bw() + facet_wrap(~ dataset) + 
  xlab("") + ylab("Longest palindrome arm length/read length") + 
  coord_cartesian(ylim = c(0, max(df$relMaxArmLengthPalindrome) + 0.005)) + 
  stat_summary(data = df %>%
                 dplyr::group_by(dataset, aligned) %>%
                 dplyr::summarize(relMaxArmLengthPalindrome = 
                                    length(relMaxArmLengthPalindrome)), 
               fun.data = function(x) {return(c(y = max(df$relMaxArmLengthPalindrome),
                                                label = x))}, 
               geom = "text", alpha = 1, size = 2.5, vjust = -1)
dev.off()

g1 <- ggplot(df, aes(x = aligned, y = maxArmLengthPalindrome)) + 
  geom_violin() + 
  theme_bw() + facet_wrap(~ dataset, nrow = 1) + 
  xlab("") + ylab("Longest palindrome arm length") + 
  coord_cartesian(ylim = c(0, max(df$maxArmLengthPalindrome) * 1.04)) + 
  stat_summary(data = df %>%
                 dplyr::group_by(dataset, aligned) %>%
                 dplyr::summarize(maxArmLengthPalindrome = 
                                    length(maxArmLengthPalindrome)), 
               fun.data = function(x) {return(c(y = max(df$maxArmLengthPalindrome),
                                                label = x))}, 
               geom = "text", alpha = 1, size = 2.5, vjust = -1)

df$hasSuppAlignment <- df$nbrSupplementaryAlignments > 0
df0 <- df %>% dplyr::filter(aligned == "Aligned")
g2 <- ggplot(df0, 
             aes(x = hasSuppAlignment, y = maxArmLengthPalindrome)) + 
  geom_violin() + 
  theme_bw() + facet_wrap(~ dataset, nrow = 1) + 
  xlab("Read has supplementary alignment") + 
  ylab("Longest palindrome arm length") + 
  coord_cartesian(ylim = c(0, max(df0$maxArmLengthPalindrome) * 1.04)) + 
  stat_summary(data = df0 %>%
                 dplyr::group_by(dataset, hasSuppAlignment) %>%
                 dplyr::summarize(maxArmLengthPalindrome = 
                                    length(maxArmLengthPalindrome)), 
               fun.data = function(x) {return(c(y = max(df0$maxArmLengthPalindrome),
                                                label = x))}, 
               geom = "text", alpha = 1, size = 2.5, vjust = -1)

df0 <- df %>% dplyr::filter(aligned == "Aligned") %>%
  dplyr::filter(dataset %in% c("ONT-DCS108-HAP", "ONT-DCS108-NSK007"))
g3 <- ggplot(df0,
             aes(x = strands_dist0, y = maxArmLengthPalindrome)) + 
    geom_violin(aes(fill = strands_dist0), alpha = 0.5) + 
    theme_bw() + facet_wrap(~ dataset) + 
    xlab("") + ylab("Longest palindrome arm length") + 
    scale_fill_manual(values = catcols, name = "") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "none") + 
    coord_cartesian(ylim = c(0, max(df0$maxArmLengthPalindrome) * 1.04)) + 
    stat_summary(data = df0 %>%
                   dplyr::group_by(dataset, strands_dist0) %>%
                   dplyr::summarize(maxArmLengthPalindrome = 
                                      length(maxArmLengthPalindrome)), 
                 fun.data = function(x) {return(c(y = max(df0$maxArmLengthPalindrome),
                                                  label = x))}, 
                 geom = "text", alpha = 1, size = 2.5, vjust = -1)

g4 <- ggplot(df0,
             aes(x = strands_dist0, y = relMaxArmLengthPalindrome)) + 
    geom_violin(aes(fill = strands_dist0), alpha = 0.5) + 
    theme_bw() + facet_wrap(~ dataset) + 
    xlab("") + ylab("Longest palindrome arm length/read length") + 
    scale_fill_manual(values = catcols, name = "") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "none") + 
    coord_cartesian(ylim = c(0, max(df0$relMaxArmLengthPalindrome) * 1.04)) + 
    stat_summary(data = df0 %>%
                   dplyr::group_by(dataset, strands_dist0) %>%
                   dplyr::summarize(relMaxArmLengthPalindrome = 
                                      length(relMaxArmLengthPalindrome)), 
                 fun.data = function(x) {return(c(y = max(df0$relMaxArmLengthPalindrome),
                                                  label = x))}, 
                 geom = "text", alpha = 1, size = 2.5, vjust = -1)

png(gsub("\\.rds$", "_summary.png", outrds),
    width = 8, height = 13, unit = "in", res = 400)
cowplot::plot_grid(
  g1, 
  g2,
  cowplot::plot_grid(g3, g4, nrow = 1, labels = c("C", "D"),
                     align = "h", axis = "bt"),
  ncol = 1, labels = c("A", "B", ""), 
  rel_heights = c(1, 1, 1.7)
)
dev.off()

png(gsub("\\.rds$", "_summary_sub.png", outrds),
    width = 8, height = 6, unit = "in", res = 400)
cowplot::plot_grid(
  g1, 
  g2,
  ncol = 1, labels = c("A", "B"), 
  rel_heights = c(1, 1)
)
dev.off()

png(gsub("\\.rds$", "_summary_sub2.png", outrds),
    width = 8, height = 3, unit = "in", res = 400)
print(g2)
dev.off()

df$spacing <- df$totalPalindromeLength - 2*df$maxArmLengthPalindrome
df0 <- df %>% dplyr::filter(aligned == "Aligned") %>%
  dplyr::filter(dataset %in% c("ONT-DCS108-HAP", "ONT-DCS108-NSK007"))
png(gsub("\\.rds$", "_stranddist_vs_spacing.png", outrds),
    width = 4, height = 7, unit = "in", res = 400)
ggplot(df0,
       aes(x = strands_dist0, y = spacing)) + 
  geom_violin(aes(fill = strands_dist0), alpha = 0.5) + 
  theme_bw() + facet_wrap(~ dataset) + 
  xlab("") + ylab("Distance between palindrome arms") + 
  scale_fill_manual(values = catcols, name = "") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none") + 
  coord_cartesian(ylim = c(0, max(df0$spacing) + 50)) + 
  stat_summary(data = df0 %>%
                 dplyr::group_by(dataset, strands_dist0) %>%
                 dplyr::summarize(spacing = length(spacing)), 
               fun.data = function(x) {return(c(y = max(df0$spacing),
                                                label = x))}, 
               geom = "text", alpha = 1, size = 2.5, vjust = -1)
dev.off()

saveRDS(df, file = outrds)
date()
sessionInfo()


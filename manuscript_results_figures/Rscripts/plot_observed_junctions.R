args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2, lib.loc = "/home/charlotte/R/x86_64-pc-linux-gnu-library/3.5")
  library(GenomicRanges)
  library(Biostrings)
  library(BSgenome)
  library(BSgenome.Hsapiens.Ensembl.GRCh3890primary, lib.loc = "/home/Shared/data/seq/hussain_bath_nanopore_rnaseq/reference/BSgenome/rlib")
})
source("manuscript_results_figures/Rscripts/remap_sample_names.R")

datasets <- strsplit(datasets, ",")[[1]]
datasets <- c(datasets, "Illumina")
names(datasets) <- datasets
conditions <- strsplit(conditions, ",")[[1]]

print(datasets)
print(conditions)
print(ntxeq)
print(outrds)

genome <- BSgenome.Hsapiens.Ensembl.GRCh3890primary

## Determine the "canonical" donor/acceptor pairs
dn <- c("AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT",
        "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT")
allDAs <- unlist(lapply(dn, function(i) paste0(i, dn)))
revcompwhere <- sapply(allDAs, function(i) {
  which(allDAs == as.character(reverseComplement(DNAString(i))))
})
canonicalDAs <- allDAs[revcompwhere <= seq_along(revcompwhere)]

harmonizeDonorAcceptor <- function(ddaa, canonicalDAs) {
  if (ddaa %in% canonicalDAs) {
    ddaa
  } else {
    tmp <- as.character(reverseComplement(DNAString(ddaa)))
    if (tmp %in% canonicalDAs) {
      tmp
    } else {
      "NNNN"
    }
  }
}

illumina_files <- list.files("Illumina/STAR", 
                             pattern = "_Aligned.sortedByCoord.out_junctions.rds", 
                             recursive = TRUE, full.names = TRUE)
illumina_files <- illumina_files[gsub("_Aligned.sortedByCoord.out_junctions.rds", "", 
                                      basename(illumina_files)) %in% 
                                   sample_annotation$sample_orig[sample_annotation$condition %in% 
                                                                   conditions]]
all_illumina_juncs <- do.call(c, lapply(illumina_files, function(f) {
  readRDS(f)
}))

juncs <- do.call(dplyr::bind_rows, lapply(datasets, function(ds) {
  print(paste0(ds, ":"))
  if (ds == "Illumina") {
    files <- list.files(paste0(ds, "/STAR"), 
                        pattern = "_Aligned.sortedByCoord.out_junctions.rds", 
                        recursive = TRUE, full.names = TRUE)
    names(files) <- gsub("_Aligned.sortedByCoord.out_junctions.rds", 
                         "", basename(files))
    files <- files[names(files) %in% 
                     sample_annotation$sample_orig[sample_annotation$condition %in% conditions]]
  } else {
    files <- list.files(paste0(ds, "/minimap2genome"), 
                        pattern = "_minimap_genome_s_junctions.rds", 
                        recursive = TRUE, full.names = TRUE)
    names(files) <- gsub("_minimap_genome_s_junctions.rds", 
                         "", basename(files))
    files <- files[names(files) %in% 
                     sample_annotation$sample_orig[sample_annotation$condition %in% conditions]]
  }
  do.call(dplyr::bind_rows, lapply(names(files), function(nm) {
    print(paste0("   ", nm))
    f <- readRDS(files[nm])
    m <- intersect(seqlevels(f), seqlevels(genome))
    f <- subset(f, seqnames %in% m)
    ## Remove introns shorter than 4 nt
    f <- f[width(f) >= 4]
    gs <- suppressWarnings(getSeq(genome, f))
    donors <- as.character(Biostrings::subseq(gs, 1, 2))
    acceptors <- as.character(Biostrings::subseq(gs, width(gs) - 1, width(gs)))
    data.frame(nbrCovReads = mcols(f)$score,
               distToClosestRefJunc = mcols(f)$distToClosestRefJunc,
               detByIllumina = f %in% all_illumina_juncs, 
               donorAcceptor = vapply(paste0(donors, acceptors), 
                                      function(w) {
                                        harmonizeDonorAcceptor(w, canonicalDAs)
                                      },
                                      "NNNN"),
               sample = remap[nm],
               dataset = remapds[ds],
               stringsAsFactors = FALSE)
  }))
}))

convertJType <- c(annotatedExact = "Identical to annotated junction",
                  annotatedLessThan5 = "Not identical, but less than 5bp from annotated junction",
                  annotatedMoreThan5 = "More than 5bp from annotated junction")

jTypeLevels <- c("More than 5bp from annotated junction",
                 "Not identical, but less than 5bp from annotated junction",
                 "Identical to annotated junction")

jTypeCols <- c(`Identical to annotated junction` = "#90C987",
               `Not identical, but less than 5bp from annotated junction` = "#B17BA6",
               `More than 5bp from annotated junction` = "#7BAFDE")

## -------------------------------------------------------------------------- ##
## Donor/acceptor pairs
## -------------------------------------------------------------------------- ##
getDonorAcceptorDf <- function(juncs, minCoverage) {
  juncs %>% dplyr::filter(nbrCovReads >= minCoverage) %>%
    dplyr::mutate(jType = ifelse(distToClosestRefJunc == 0, "annotatedExact",
                                 ifelse(distToClosestRefJunc <= 5, 
                                        "annotatedLessThan5", "annotatedMoreThan5"))) %>%
    dplyr::group_by(sample, dataset, jType, donorAcceptor) %>%
    dplyr::tally() %>%
    dplyr::mutate(donorAcceptor = factor(donorAcceptor, levels = c(setdiff(unique(juncs$donorAcceptor), "GTAG"), "GTAG"))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(jType = convertJType[jType]) %>%
    dplyr::mutate(jType = factor(jType, levels = jTypeLevels)) %>%
    dplyr::mutate(dataset = factor(dataset, levels = c(setdiff(levels(factor(dataset)), 
                                                               "Illumina"), "Illumina")))
}

plots <- list()

for (minCov in c(1, 5)) {
  df0 <- getDonorAcceptorDf(juncs, minCov)
  
  gglayers <- list(
    geom_bar(stat = "identity", position = "fill", alpha = 0.5),
    facet_grid(jType ~ dataset, scales = "free_x", space = "free_x"),
    theme_bw(),
    xlab(""),
    theme(strip.text = element_text(size = 8)),
    scale_y_continuous(expand = c(0, 0, 0.05, 0)),
    ylab("Fraction of junctions with indicated donor/acceptor sequence"),
    ggtitle(paste0("Junctions supported by at least ", minCov, 
                   ifelse(minCov == 1, " read", " reads")))
  )
  
  plots[[paste0("donoracceptor_minCov", minCov)]] <- 
    ggplot(df0 %>% dplyr::mutate(sample = removeDatasetFromSample(sample, dataset)),
           aes(x = sample, y = n, fill = donorAcceptor)) + 
    gglayers
  
  plots[[paste0("donoracceptor_twogroups_minCov", minCov)]] <- 
    ggplot(df0 %>% dplyr::mutate(sample = removeDatasetFromSample(sample, dataset)) %>% 
             dplyr::mutate(donorAcceptor = as.character(donorAcceptor)) %>%
             dplyr::mutate(donorAcceptor = replace(donorAcceptor, 
                                                   donorAcceptor != "GTAG", "other")) %>%
             dplyr::mutate(donorAcceptor = factor(donorAcceptor, 
                                                  levels = c("other", "GTAG"))),
           aes(x = sample, y = n, fill = donorAcceptor)) + 
    gglayers + 
    scale_fill_manual(values = c(GTAG = "red", other = "darkgreen"), 
                      name = "Donor/acceptor site") + 
    theme(legend.position = "bottom")
  
  png(gsub("\\.rds$", paste0("junctions_with_atleast_", minCov, 
                             "_reads_donoracceptor.png"), outrds), 
      width = 24, height = 15, unit = "in", res = 400)
  print(plots[[paste0("donoracceptor_minCov", minCov)]])
  dev.off()
  
  png(gsub("\\.rds$", paste0("junctions_with_atleast_", minCov, 
                             "_reads_donoracceptor_twogroups.png"), outrds), 
      width = 24, height = 15, unit = "in", res = 400)
  print(plots[[paste0("donoracceptor_twogroups_minCov", minCov)]])
  dev.off()
}

## -------------------------------------------------------------------------- ##
## Agreement with annotation
## -------------------------------------------------------------------------- ##
getJuncDf <- function(juncs, minCoverage) {
  juncs %>% dplyr::filter(nbrCovReads >= minCoverage) %>%
    dplyr::group_by(dataset, sample) %>%
    dplyr::summarize(annotatedExact = sum(distToClosestRefJunc == 0),
                     annotatedLessThan5 = sum(distToClosestRefJunc <= 5 & 
                                                distToClosestRefJunc > 0),
                     annotatedMoreThan5 = sum(distToClosestRefJunc > 5)) %>%
    tidyr::gather(jType, nbrJunctions, -dataset, -sample) %>%
    dplyr::mutate(jType = convertJType[jType]) %>%
    dplyr::mutate(jType = factor(jType, levels = jTypeLevels)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(dataset = factor(dataset, levels = c(setdiff(levels(factor(dataset)), 
                                                               "Illumina"), "Illumina")))
}

for (minCov in c(1, 5)) {
  df <- getJuncDf(juncs, minCov)
  
  plots[[paste0("junctions_minCov", minCov)]] <- 
    ggplot(df %>% dplyr::mutate(sample = removeDatasetFromSample(sample, dataset)),
           aes(x = sample, y = nbrJunctions, fill = jType, color = jType)) + 
    geom_bar(stat = "identity", position = "fill") + 
    theme_bw() + 
    scale_y_continuous(expand = c(0, 0, 0.05, 0)) + 
    facet_grid(~ dataset, scales = "free_x", space = "free_x") + 
    theme(legend.position = "bottom",
          strip.text = element_text(size = 8)) + 
    xlab("") + ylab("Fraction of junctions") + 
    scale_fill_manual(name = "", values = jTypeCols) + 
    scale_color_manual(name = "", values = jTypeCols) + 
    ggtitle(paste0("Junctions supported by at least ", minCov, 
                   ifelse(minCov == 1, " read", " reads")))
  
  png(gsub("\\.rds$", paste0("junctions_with_atleast_", minCov, "reads.png"), 
           outrds), width = 16, height = 5.5, 
      unit = "in", res = 400)
  print(plots[[paste0("junctions_minCov", minCov)]])
  dev.off()
}

## -------------------------------------------------------------------------- ##
## Observed in Illumina
## -------------------------------------------------------------------------- ##
getIlmnJuncDf <- function(juncs, minCoverage) {
  juncs %>% dplyr::filter(nbrCovReads >= minCoverage) %>%
    dplyr::group_by(dataset, sample, detByIllumina) %>%
    dplyr::summarize(annotatedExact = sum(distToClosestRefJunc == 0),
                     annotatedLessThan5 = sum(distToClosestRefJunc <= 5 & 
                                                distToClosestRefJunc > 0),
                     annotatedMoreThan5 = sum(distToClosestRefJunc > 5)) %>%
    tidyr::gather(jType, nbrJunctions, -dataset, -sample, -detByIllumina) %>%
    dplyr::mutate(jType = convertJType[jType]) %>%
    dplyr::mutate(jType = factor(jType, levels = jTypeLevels)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(dataset = factor(dataset, levels = c(setdiff(levels(factor(dataset)), 
                                                               "Illumina"), "Illumina")))
}

for (minCov in c(1, 5)) {
  dfilmn <- getIlmnJuncDf(juncs, minCov)
  
  plots[[paste0("detbyillumina_minCov", minCov)]] <- 
    ggplot(dfilmn %>% 
             dplyr::mutate(sample = removeDatasetFromSample(sample, dataset)) %>%
             dplyr::filter(dataset != "Illumina"), 
           aes(x = sample, y = nbrJunctions, fill = detByIllumina, 
               color = detByIllumina)) + 
    geom_bar(stat = "identity", position = "fill", alpha = 0.5) + 
    theme_bw() + 
    scale_y_continuous(expand = c(0, 0, 0.05, 0)) + 
    facet_grid(jType ~ dataset, scales = "free_x", space = "free_x") + 
    theme(legend.position = "bottom",
          strip.text = element_text(size = 8)) + 
    xlab("") + ylab("Fraction of junctions") + 
    scale_fill_manual(
      name = "Detected in any Illumina sample", 
      values = c(`TRUE` = "darkgreen", `FALSE` = "red")) + 
    scale_color_manual(
      name = "Detected in any Illumina sample", 
      values = c(`TRUE` = "darkgreen", `FALSE` = "red")) + 
    ggtitle(paste0("Junctions supported by at least ", minCov, 
                   ifelse(minCov == 1, " read", " reads")))
  
  png(gsub("\\.rds$", paste0("_junctions_with_atleast_", minCov, "_reads_detbyillumina.png"),
           outrds), width = 24, height = 15, 
      unit = "in", res = 400)
  print(plots[[paste0("detbyillumina_minCov", minCov)]])
  dev.off()
}

## -------------------------------------------------------------------------- ##
## Combined plots
## -------------------------------------------------------------------------- ##
png(gsub("\\.rds$", "_detbyillumina_donoracceptor_mincov5.png", outrds), 
    width = 12, height = 20, unit = "in", res = 400)
print(cowplot::plot_grid(
  plots[["detbyillumina_minCov5"]],
  plots[["donoracceptor_twogroups_minCov5"]],
  ncol = 1, labels = c("A", "B"), rel_heights = c(1, 1)
))
dev.off()

ntxeq <- readRDS(ntxeq)

png(gsub("\\.rds$", "_junctions_minCov1_5.png", outrds),
    width = 12, height = 12, unit = "in", res = 400)
print(cowplot::plot_grid(
  plots[["junctions_minCov1"]] + theme(legend.position = "none") + 
    ylab("Fraction of\njunctions"),
  plots[["junctions_minCov5"]] + theme(legend.position = "none") + 
    ylab("Fraction of\njunctions"),
  cowplot::get_legend(plots[["junctions_minCov5"]]),
  ntxeq$ptxeqfull + ylab("Number of transcripts\nin equivalence class"),
  ntxeq$ptxeqzoom + ylab("Number of transcripts\nin equivalence class"), 
  ncol = 1, labels = c("A", "B", "", "C", "D"), rel_heights = c(1, 1, 0.2, 1, 1)
))

saveRDS(NULL, file = outrds)
date()
sessionInfo()

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2, lib.loc = "/home/charlotte/R/x86_64-pc-linux-gnu-library/3.5")
  library(cowplot)
  library(stringr)
})
source("manuscript_results_figures/Rscripts/remap_sample_names.R")

datasets <- strsplit(datasets, ",")[[1]]
names(datasets) <- datasets
conditions <- strsplit(conditions, ",")[[1]]

print(datasets)
print(conditions)
print(outrds)

## Read info
readInfo <- lapply(datasets, function(ds) {
  rtbl <- read.delim(paste0(ds, "/output/", ds, "_nbr_reads.txt"), 
                     header = TRUE, as.is = TRUE)
  rtbl[!grepl("_orig", rtbl$sample), ] %>%
    dplyr::mutate(dataset = remapds[ds]) %>%
    dplyr::mutate(sample = remap[sample])
})

## Number of reads
nReads <- do.call(dplyr::bind_rows, readInfo) %>%
  dplyr::left_join(sample_annotation %>% dplyr::select(sample_remap, condition),
                   by = c("sample" = "sample_remap")) %>%
  dplyr::filter(condition %in% conditions) %>%
  dplyr::select(-condition)

## Abundance info
abundanceInfo <- lapply(datasets, function(ds) {
  readRDS(paste0(ds, "/output/", ds, "_all_abundances.rds"))
})

## Keep only the desired abundance measures
abundanceInfo <- lapply(abundanceInfo, function(w) {
  w$gene_abundances <- 
    w$gene_abundances[, grep("count__salmon31|count__salmonminimap2_p0.99|count__wubminimap2|count__featurecountsminimap2primary", 
                             colnames(w$gene_abundances))]
  w
})

genelevel <- Reduce(function(...) dplyr::full_join(..., by = "gene"), 
                    lapply(abundanceInfo, function(ab) {
                      ab$gene_abundances %>% tibble::rownames_to_column("gene")
                    })) %>%
  tibble::column_to_rownames("gene")
totCounts <- 
  data.frame(totCount = colSums(genelevel), method = colnames(genelevel), 
             stringsAsFactors = FALSE) %>%
  tidyr::separate(method, into = c("sample", "type", "method"), sep = "__") %>%
  dplyr::mutate(sample = remap[sample]) %>% 
  dplyr::left_join(sample_annotation %>% dplyr::select(sample_remap, condition), 
                   by = c("sample" = "sample_remap")) %>%
  dplyr::filter(condition %in% conditions) %>% 
  dplyr::select(sample, method, totCount) %>%
  tidyr::spread(method, totCount)

colvec <- c(`Total number of reads` = "#777777",
            `Number of reads with a primary alignment to the genome` = "#E8601C",
            `Number of reads with a primary alignment to the transcriptome` = "#90C987",
            `Number of assigned reads, salmonminimap2_p0.99` = "#1965B0",
            `Number of assigned reads, fCminimap2primary` = "#882E72",
            `Number of assigned reads, salmon31` = "#DC050C",
            `Number of assigned reads, wubminimap2` = "#55A1B1")

nReadsData <- nReads %>% 
  dplyr::select(sample, dataset, nReads, nAlignedReadsGenome,
                nAlignedReadsTxome) %>%
  dplyr::full_join(totCounts, by = "sample") %>%
  tidyr::gather(rtype, nReads, -sample, -dataset) %>%
  dplyr::mutate(
    rtype = replace(rtype, rtype == "nReads",
                    "Total number of reads"),
    rtype = replace(rtype, rtype == "nAlignedReadsGenome",
                    "Number of reads with a primary alignment to the genome"),
    rtype = replace(rtype, rtype == "nAlignedReadsTxome",
                    "Number of reads with a primary alignment to the transcriptome"),
    rtype = replace(rtype, rtype == "featurecountsminimap2primary",
                    "Number of assigned reads, fCminimap2primary"),
    rtype = replace(rtype, rtype == "salmon31",
                    "Number of assigned reads, salmon31"),
    rtype = replace(rtype, rtype == "salmonminimap2_p0.99",
                    "Number of assigned reads, salmonminimap2_p0.99"),
    rtype = replace(rtype, rtype == "wubminimap2",
                    "Number of assigned reads, wubminimap2")
  ) %>%
  dplyr::mutate(
    rtype = factor(rtype, levels = names(colvec)
    ))

png(gsub("\\.rds$", ".png", outrds), width = 16, height = 6, unit = "in", res = 400)
print(ggplot(nReadsData %>% 
               dplyr::mutate(nReads = nReads/1e6) %>%
               dplyr::mutate(sample = removeDatasetFromSample(sample, dataset)),
             aes(x = sample, y = nReads, fill = rtype)) + 
        geom_bar(stat = "identity", position = "dodge") + theme_bw() + 
        scale_y_continuous(expand = c(0, 0, 0.05, 0), 
                           limits = c(0, max(nReadsData$nReads)/1e6)) + 
        facet_grid(~ dataset, scales = "free_x", space = "free_x") + 
        scale_fill_manual(values = colvec, name = "") +
        theme(legend.position = "bottom",
              strip.text = element_text(size = 8)) + 
        guides(fill = guide_legend(nrow = 3, byrow = TRUE)) + 
        xlab("") + ylab("Number of reads (Mio.)") + 
        stat_summary(data = nReadsData %>% 
                       dplyr::mutate(nReads = nReads/1e6) %>%
                       dplyr::mutate(sample = removeDatasetFromSample(sample, dataset)) %>%
                       dplyr::group_by(sample, dataset) %>% 
                       dplyr::mutate(
                         nReads = round(nReads/nReads[rtype=="Total number of reads"]*100)), 
                     fun.data = function(x) {
                       return(c(y = ifelse(x == 100, -100000, 1e-4), label = x))}, 
                     geom = "text", alpha = 1, size = 2, vjust = -1, 
                     position = position_dodge(width = 0.9)))

dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()

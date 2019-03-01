args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(bamsignals)
  library(dplyr)
  library(ggplot2)
  library(tximport)
  library(tibble)
})
source("manuscript_results_figures/Rscripts/remap_sample_names.R")

datasets <- strsplit(datasets, ",")[[1]]
datasets <- c(datasets, "Illumina")
conditions <- strsplit(conditions, ",")[[1]]

print(datasets)
print(conditions)
print(tx2gene)
print(gtf)
print(outrds)

tx2gene <- readRDS(tx2gene)
genes <- rtracklayer::import(gtf)
genes <- subset(genes, type == "exon")

## Get coverages for all bases in all transcripts, for all samples
baseCoverages <- do.call(dplyr::bind_rows, lapply(datasets, function(ds) {
  print(paste0(ds, ":"))
  if (ds == "Illumina") {
    files <- list.files(paste0(ds, "/STAR"), 
                        pattern = "_Aligned.sortedByCoord.out.bam$", 
                        recursive = TRUE, full.names = TRUE)
    names(files) <- gsub("_Aligned.sortedByCoord.out.bam", 
                         "", basename(files))
    files <- files[names(files) %in% 
                     sample_annotation$sample_orig[sample_annotation$condition %in% conditions]]
  } else {
    files <- list.files(paste0(ds, "/minimap2genome"), 
                        pattern = "_minimap_genome_s.bam$", 
                        recursive = TRUE, full.names = TRUE)
    names(files) <- gsub("_minimap_genome_s.bam", 
                         "", basename(files))
    files <- files[names(files) %in% 
                     sample_annotation$sample_orig[sample_annotation$condition %in% conditions]]
  }
  do.call(dplyr::bind_rows, lapply(names(files), function(nm) {
    print(paste0("   ", nm))
    sigs <- bamCoverage(files[nm], genes, verbose = FALSE)
    as.data.frame(genes) %>% 
      dplyr::mutate(totCov = sapply(sigs, sum),
                    totNBases = sapply(sigs, length),
                    totNonzero = sapply(sigs, function(w) sum(w > 0))) %>%
      dplyr::group_by(transcript_id) %>%
      dplyr::summarize(txLength = sum(width),
                       totCov = sum(totCov), 
                       totNBases = sum(totNBases),
                       totNonzero = sum(totNonzero),
                       txBiotype = unique(transcript_biotype),
                       gene = unique(gene_id),
                       geneName = unique(gene_name),
                       geneBiotype = unique(gene_biotype)) %>%
      dplyr::mutate(dataset = remapds[ds],
                    sample = remap[nm])
  }))
}))

## Add expression (TPM) in Illumina samples
files <- list.files("Illumina/salmon31", pattern = "quant.sf", 
                    recursive = TRUE, full.names = TRUE)
names(files) <- basename(dirname(files))
files <- files[names(files) %in% 
                 sample_annotation$sample_orig[sample_annotation$condition %in% conditions]]
salmon <- tximport(files = files, type = "salmon", txOut = TRUE)

tpm <- as.data.frame(salmon$abundance) %>%
  tibble::rownames_to_column("transcript_id") %>%
  dplyr::mutate(transcript_id = gsub("\\.[0-9]+$", "", transcript_id)) %>%
  tidyr::gather(key = sample, value = TPM, -transcript_id) %>%
  dplyr::mutate(sample = remap[sample],
                dataset = "Illumina")

df <- baseCoverages %>% 
  dplyr::left_join(tpm, by = c("sample", "dataset", "transcript_id"))

## Plot transcript length vs number of covered bases for Illumina data,
## at different expression cutoffs
plots <- list()
png(gsub("\\.rds$", paste0("_txlength_vs_numcovbases_illumina_minTPM1.png"), outrds),
    width = 7, height = 7, unit = "in", res = 400)
df0 <- df %>% dplyr::filter(dataset == "Illumina") %>%
  dplyr::group_by(transcript_id) %>%
  dplyr::summarize(aveTPM = mean(TPM),
                   aveNbrNonzero = mean(totNonzero),
                   txLength = mean(txLength)) %>%
  dplyr::filter(aveTPM > 1)
plots[[paste0("txlength_vs_basecoverage_minTPM1")]] <- 
  ggplot(df0, aes(x = txLength, y = aveNbrNonzero)) + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_hex(bins = 100, aes(fill = stat(density))) + 
  scale_fill_gradient(name = "", low = "bisque2", high = "darkblue") + 
  theme_bw() + 
  xlab("Transcript length") + 
  ylab("Average number of covered bases in Illumina samples") + 
  ggtitle(paste0("Transcripts with average TPM > 1"))
print(plots[[paste0("txlength_vs_basecoverage_minTPM1")]])
dev.off()

saveRDS(plots, file = outrds)
date()
sessionInfo()


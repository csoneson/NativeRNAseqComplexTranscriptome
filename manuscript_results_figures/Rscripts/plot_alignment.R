args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(stringr)
  library(GenomicAlignments)
})
source("manuscript_results_figures/Rscripts/remap_sample_names.R")

datasets <- strsplit(datasets, ",")[[1]]
names(datasets) <- datasets
conditions <- strsplit(conditions, ",")[[1]]

print(datasets)
print(conditions)
print(tx2gene)
print(duptx)  ## List of duplicated transcripts (from the Salmon index building)
print(ilmnnjunc)  ## Number of junctions covered by each Illumina read
print(primsupp)  ## Plot, primary/supplementary alignment overlap
print(txlvsbc)  ## Plot, transcript length vs base coverage (Illumina)
print(outrds)

tx2gene <- readRDS(tx2gene)

## Read info
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
  
  rdt <- rd$txomebams
  rdt <- rdt[!grepl("_orig", names(rdt))]
  rdt <- rdt[names(rdt) %in% 
               sample_annotation$sample_orig[sample_annotation$condition %in% conditions]]
  rdt <- do.call(dplyr::bind_rows, lapply(names(rdt), function(nm) {
    rdt[[nm]]$allAlignments %>% 
      dplyr::mutate(sample = remap[nm], dataset = remapds[ds], 
                    rtype = "Reads aligning to the transcriptome") %>%
      dplyr::mutate(fracM = nbrM/readLength,
                    fracS = nbrS/readLength,
                    fracI = nbrI/readLength,
                    accuracy = (nbrM - NM + nbrI + nbrD)/(nbrM + nbrI + nbrD))
  })) %>%
    dplyr::left_join(tx2gene %>% dplyr::rename(rname = tx, gname = gene) %>%
                       dplyr::select(rname, gname),
                     by = "rname")

  rdt_p0.99 <- rd$txomebams_p0.99
  rdt_p0.99 <- rdt_p0.99[!grepl("_orig", names(rdt_p0.99))]
  rdt_p0.99 <- rdt_p0.99[names(rdt_p0.99) %in% 
                           sample_annotation$sample_orig[sample_annotation$condition %in% 
                                                           conditions]]
  rdt_p0.99 <- do.call(dplyr::bind_rows, lapply(names(rdt_p0.99), function(nm) {
    rdt_p0.99[[nm]]$allAlignments %>% 
      dplyr::mutate(sample = remap[nm], dataset = remapds[ds], 
                    rtype = "Reads aligning to the transcriptome (-p 0.99)") %>%
      dplyr::mutate(fracM = nbrM/readLength,
                    fracS = nbrS/readLength,
                    fracI = nbrI/readLength,
                    accuracy = (nbrM - NM + nbrI + nbrD)/(nbrM + nbrI + nbrD))
  })) %>%
    dplyr::left_join(tx2gene %>% dplyr::rename(rname = tx, gname = gene) %>%
                       dplyr::select(rname, gname),
                     by = "rname")
  
  rtbl <- rd$nReadTable
  rtbl <- rtbl[!grepl("_orig", rtbl$sample), ] %>%
    dplyr::mutate(dataset = remapds[ds],
                  sample = remap[sample]) %>%
    dplyr::left_join(sample_annotation %>% dplyr::select(sample_remap, condition),
                     by = c("sample" = "sample_remap")) %>%
    dplyr::filter(condition %in% conditions) %>%
    dplyr::select(-condition)
  
  list(rdf = rdf, rdg = rdg, rdt = rdt, rdt_p0.99 = rdt_p0.99, rtbl = rtbl)
})

## Number of reads
nReads <- do.call(dplyr::bind_rows, lapply(readInfo, function(x) x$rtbl))

## All reads
allReads <- do.call(dplyr::bind_rows, lapply(readInfo, function(x) x$rdf))

## Genome alignments
gAlign <- do.call(dplyr::bind_rows, lapply(readInfo, function(x) x$rdg))

## Transcriptome alignments
tAlign <- do.call(dplyr::bind_rows, lapply(readInfo, function(x) x$rdt))

## Transcriptome alignments, -p 0.99
tAlign_p0.99 <- do.call(dplyr::bind_rows, lapply(readInfo, function(x) x$rdt_p0.99))

## Genome alignments for reads that do not align to the transcriptome
gnottAlign <- gAlign %>%
  dplyr::filter(!(read %in% tAlign_p0.99$read)) %>%
  dplyr::mutate(rtype = "Reads aligning to the genome but not the transcriptome")

## Help function to calculate the N50 of a vector of values
calcN50 <- function(x) {
  ## x = vector with read/alignment/contig lengths
  x <- sort(x, decreasing = TRUE)
  totLength <- sum(x)
  data.frame(lengths = x, cumLengths = cumsum(x)) %>%
    dplyr::filter(cumLengths > totLength/2) %>%
    dplyr::filter(cumLengths == min(cumLengths)) %>%
    dplyr::pull(lengths)
}

## Calculate summary statistics for primary genome alignments
dfsumm <- dplyr::bind_rows(
  gAlign %>% dplyr::filter(flag %in% c(0, 16)) %>%
    dplyr::group_by(sample, dataset) %>%
    dplyr::summarize(totReadLengthAlignedGenomePrimary = sum(readLength),
                     totAlignedLengthGenomePrimary = sum(alignedLength),
                     medianReadLengthAlignedGenomePrimary = median(readLength),
                     medianAlignedLengthGenomePrimary = median(alignedLength),
                     maxAlignedLengthGenomePrimary = max(alignedLength),
                     N50ReadLengthAlignedGenomePrimary = calcN50(readLength),
                     N50AlignedLengthGenomePrimary = calcN50(alignedLength)),
  gAlign %>% dplyr::filter(flag %in% c(0, 16)) %>%
    dplyr::summarize(totReadLengthAlignedGenomePrimary = sum(readLength),
                     totAlignedLengthGenomePrimary = sum(alignedLength),
                     medianReadLengthAlignedGenomePrimary = median(readLength),
                     medianAlignedLengthGenomePrimary = median(alignedLength),
                     maxAlignedLengthGenomePrimary = max(alignedLength),
                     N50ReadLengthAlignedGenomePrimary = calcN50(readLength),
                     N50AlignedLengthGenomePrimary = calcN50(alignedLength)) %>%
    dplyr::mutate(sample = "Overall", dataset = "Overall"),
  gAlign %>% dplyr::filter(flag %in% c(0, 16) & 
                             dataset %in% c("ONT-RNA001-HAP", "ONT-RNA001-HEK")) %>%
    dplyr::summarize(totReadLengthAlignedGenomePrimary = sum(readLength),
                     totAlignedLengthGenomePrimary = sum(alignedLength),
                     medianReadLengthAlignedGenomePrimary = median(readLength),
                     medianAlignedLengthGenomePrimary = median(alignedLength),
                     maxAlignedLengthGenomePrimary = max(alignedLength),
                     N50ReadLengthAlignedGenomePrimary = calcN50(readLength),
                     N50AlignedLengthGenomePrimary = calcN50(alignedLength)) %>%
    dplyr::mutate(sample = "nativeRNA", dataset = "nativeRNA"),
  gAlign %>% dplyr::filter(flag %in% c(0, 16)) %>%
    dplyr::group_by(dataset) %>%
    dplyr::summarize(totReadLengthAlignedGenomePrimary = sum(readLength),
                     totAlignedLengthGenomePrimary = sum(alignedLength),
                     medianReadLengthAlignedGenomePrimary = median(readLength),
                     medianAlignedLengthGenomePrimary = median(alignedLength),
                     maxAlignedLengthGenomePrimary = max(alignedLength),
                     N50ReadLengthAlignedGenomePrimary = calcN50(readLength),
                     N50AlignedLengthGenomePrimary = calcN50(alignedLength)) %>%
    dplyr::mutate(sample = "AllSamples")
)

## Write alignment rates to file
write.table(dplyr::full_join(
  dplyr::bind_rows(
    nReads %>% 
      dplyr::mutate(alignmentRateGenome = nAlignedReadsGenome/nReads,
                    alignmentRateTxome = nAlignedReadsTxome/nReads,
                    alignmentRateTxome_p0.99 = nAlignedReadsTxome_p0.99/nReads) %>%
      dplyr::select(sample, dataset, nReads, nAlignedReadsGenome, 
                    alignmentRateGenome, nReadsWithSecondaryAlignmentsGenome,
                    nReadsWithSupplementaryAlignmentsGenome, 
                    nAlignedReadsTxome, alignmentRateTxome,
                    nReadsWithSecondaryAlignmentsTxome,
                    nReadsWithSupplementaryAlignmentsTxome,
                    nAlignedReadsTxome_p0.99, alignmentRateTxome_p0.99,
                    nReadsWithSecondaryAlignmentsTxome_p0.99,
                    nReadsWithSupplementaryAlignmentsTxome_p0.99),
    nReads %>% 
      dplyr::group_by(dataset) %>%
      dplyr::summarize(nReads = sum(nReads), 
                       nAlignedReadsGenome = sum(nAlignedReadsGenome),
                       nReadsWithSecondaryAlignmentsGenome = 
                         sum(nReadsWithSecondaryAlignmentsGenome),
                       nReadsWithSupplementaryAlignmentsGenome = 
                         sum(nReadsWithSupplementaryAlignmentsGenome),
                       nAlignedReadsTxome = sum(nAlignedReadsTxome),
                       nReadsWithSecondaryAlignmentsTxome = 
                         sum(nReadsWithSecondaryAlignmentsTxome),
                       nReadsWithSupplementaryAlignmentsTxome = 
                         sum(nReadsWithSupplementaryAlignmentsTxome),
                       nAlignedReadsTxome_p0.99 = sum(nAlignedReadsTxome_p0.99),
                       nReadsWithSecondaryAlignmentsTxome_p0.99 = 
                         sum(nReadsWithSecondaryAlignmentsTxome_p0.99),
                       nReadsWithSupplementaryAlignmentsTxome_p0.99 = 
                         sum(nReadsWithSupplementaryAlignmentsTxome_p0.99)) %>%
      dplyr::mutate(alignmentRateGenome = nAlignedReadsGenome/nReads,
                    alignmentRateTxome = nAlignedReadsTxome/nReads,
                    alignmentRateTxome_p0.99 = nAlignedReadsTxome_p0.99/nReads) %>%
      dplyr::mutate(sample = "AllSamples") %>%
      dplyr::select(sample, dataset, nReads, nAlignedReadsGenome, 
                    alignmentRateGenome, nReadsWithSecondaryAlignmentsGenome,
                    nReadsWithSupplementaryAlignmentsGenome, 
                    nAlignedReadsTxome, alignmentRateTxome,
                    nReadsWithSecondaryAlignmentsTxome,
                    nReadsWithSupplementaryAlignmentsTxome,
                    nAlignedReadsTxome_p0.99, alignmentRateTxome_p0.99,
                    nReadsWithSecondaryAlignmentsTxome_p0.99,
                    nReadsWithSupplementaryAlignmentsTxome_p0.99)
  ), 
  dfsumm),
  file = gsub("\\.rds$", "_alignment_rates.txt", outrds),
  sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

## ========================================================================== ##
## Number of reads with primary alignment to genome/transcriptome
## ========================================================================== ##
nReadsData <- nReads %>% 
  dplyr::select(sample, dataset, nReads, nAlignedReadsGenome,
                nAlignedReadsTxome) %>%
  tidyr::gather(rtype, nReads, -sample, -dataset) %>%
  dplyr::mutate(
    rtype = replace(rtype, rtype == "nReads",
                    "Total number of reads"),
    rtype = replace(rtype, rtype == "nAlignedReadsGenome",
                    "Number of reads with a primary alignment to the genome"),
    rtype = replace(rtype, rtype == "nAlignedReadsTxome",
                    "Number of reads with a primary alignment to the transcriptome")
  ) %>%
  dplyr::mutate(
    rtype = factor(
      rtype, 
      levels = c("Total number of reads",
                 "Number of reads with a primary alignment to the genome",
                 "Number of reads with a primary alignment to the transcriptome")
    )) %>%
  dplyr::mutate(dataset = factor(dataset, levels = ds_order[ds_order %in% dataset]))

## Same, but aggregated by dataset
nReadsDatads <- nReads %>% 
  dplyr::select(sample, dataset, nReads, nAlignedReadsGenome,
                nAlignedReadsTxome) %>%
  tidyr::gather(rtype, nReads, -sample, -dataset) %>%
  dplyr::group_by(rtype, dataset) %>%
  dplyr::summarize(nReads = sum(nReads)) %>% 
  dplyr::ungroup() %>%
  dplyr::mutate(
    rtype = replace(rtype, rtype == "nReads",
                    "Total number of reads"),
    rtype = replace(rtype, rtype == "nAlignedReadsGenome",
                    "Number of reads with a primary alignment to the genome"),
    rtype = replace(rtype, rtype == "nAlignedReadsTxome",
                    "Number of reads with a primary alignment to the transcriptome")
  ) %>%
  dplyr::mutate(
    rtype = factor(
      rtype, 
      levels = c("Total number of reads",
                 "Number of reads with a primary alignment to the genome",
                 "Number of reads with a primary alignment to the transcriptome")
    )) %>%
  dplyr::mutate(dataset = factor(dataset, levels = ds_order[ds_order %in% dataset]))

pnreads <- 
  ggplot(nReadsData %>% 
           dplyr::mutate(nReads = nReads/1e6) %>%
           dplyr::mutate(sample = removeDatasetFromSample(sample, dataset)),
         aes(x = sample, y = nReads, fill = rtype)) + 
  geom_bar(stat = "identity", position = "dodge") + theme_bw() + 
  scale_y_continuous(expand = c(0, 0, 0.05, 0)) + 
  facet_grid(~ dataset, scales = "free_x", space = "free_x") + 
  scale_fill_manual(
    values = c(
      `Total number of reads` = "#777777", 
      `Number of reads with a primary alignment to the genome` = "#E8601C", 
      `Number of reads with a primary alignment to the transcriptome` = "#90C987"), 
    name = "") + 
  theme(legend.position = "bottom",
        strip.text = element_text(size = 8)) + 
  xlab("") + ylab("Number of reads (Mio.)") + 
  stat_summary(
    data = nReadsData %>% 
      dplyr::mutate(nReads = nReads/1e6) %>%
      dplyr::mutate(sample = removeDatasetFromSample(sample, dataset)) %>%
      dplyr::group_by(sample, dataset) %>% 
      dplyr::mutate(nReads = round(nReads/nReads[rtype == "Total number of reads"]*100)), 
    fun.data = function(x) {return(c(y = ifelse(x == 100, -100000, 1e-4), label = x))}, 
    geom = "text", alpha = 1, size = 2, vjust = -1, 
    position = position_dodge(width = 0.9)) + 
  ylim(0, NA)
  
pnreadsds <- 
  ggplot(nReadsDatads %>% 
           dplyr::mutate(nReads = nReads/1e6),
         aes(x = dataset, y = nReads, fill = rtype)) + 
  geom_bar(stat = "identity", position = "dodge") + theme_bw() + 
  scale_y_continuous(expand = c(0, 0, 0.05, 0)) + 
  scale_fill_manual(
    values = c(
      `Total number of reads` = "#777777", 
      `Number of reads with a primary alignment to the genome` = "#E8601C", 
      `Number of reads with a primary alignment to the transcriptome` = "#90C987"), 
    name = "") + 
  theme(legend.position = "bottom",
        strip.text = element_text(size = 8),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  xlab("") + ylab("Number of reads (Mio.)") + 
  guides(fill = guide_legend(nrow = 3, byrow = TRUE)) + 
  stat_summary(
    data = nReadsDatads %>% 
      dplyr::mutate(nReads = nReads/1e6) %>%
      dplyr::group_by(dataset) %>% 
      dplyr::mutate(nReads = round(nReads/nReads[rtype == "Total number of reads"]*100)), 
    fun.data = function(x) {return(c(y = ifelse(x == 100, -100000, 1e-4), label = x))}, 
    geom = "text", alpha = 1, size = 4, vjust = -1, 
    position = position_dodge(width = 0.9)) + 
  ylim(0, NA)

## ========================================================================== ##
## Fraction of reads with primary alignment that also have secondary/
## supplementary alignments
## ========================================================================== ##
nReadsDataFS <- nReads %>%
  dplyr::mutate(fracSecGenome = nReadsWithSecondaryAlignmentsGenome /
                  nReadsWithPrimaryAlignmentsGenome,
                fracSuppGenome = nReadsWithSupplementaryAlignmentsGenome /
                  nReadsWithPrimaryAlignmentsGenome,
                fracSecTxome = nReadsWithSecondaryAlignmentsTxome /
                  nReadsWithPrimaryAlignmentsTxome,
                fracSuppTxome = nReadsWithSupplementaryAlignmentsTxome /
                  nReadsWithPrimaryAlignmentsTxome,
                fracSecTxome_p0.99 = nReadsWithSecondaryAlignmentsTxome_p0.99 /
                  nReadsWithPrimaryAlignmentsTxome_p0.99) %>%
  dplyr::select(sample, dataset, fracSecGenome, fracSuppGenome,
                fracSecTxome, fracSuppTxome, fracSecTxome_p0.99) %>%
  tidyr::gather(ftype, fracReads, -sample, -dataset) %>%
  dplyr::mutate(ref = stringr::str_extract(ftype, "Genome|Txome")) %>%
  dplyr::mutate(ref = replace(ref, ref == "Txome", "Transcriptome")) %>%
  dplyr::mutate(ftype = gsub("Genome|Txome", "", ftype)) %>%
  dplyr::mutate(
    ftype = replace(ftype, ftype == "fracSec", 
                    paste0("Fraction of reads with primary alignment", 
                           " that also have secondary alignment(s)")),
    ftype = replace(ftype, ftype == "fracSupp",
                    paste0("Fraction of reads with primary alignment",
                           " that also have supplementary alignment(s)")),
    ftype = replace(ftype, ftype == "fracSec_p0.99", 
                    paste0("Fraction of reads with primary alignment",
                           " that also have secondary alignment(s)", 
                           ", with -p=0.99"))
  )

## Summarized by dataset
nReadsDataFSds <- nReads %>%
  dplyr::group_by(dataset) %>%
  dplyr::summarize(
    fracSecGenome = sum(nReadsWithSecondaryAlignmentsGenome) /
      sum(nReadsWithPrimaryAlignmentsGenome),
    fracSuppGenome = sum(nReadsWithSupplementaryAlignmentsGenome) /
      sum(nReadsWithPrimaryAlignmentsGenome),
    fracSecTxome = sum(nReadsWithSecondaryAlignmentsTxome) /
      sum(nReadsWithPrimaryAlignmentsTxome),
    fracSuppTxome = sum(nReadsWithSupplementaryAlignmentsTxome) /
      sum(nReadsWithPrimaryAlignmentsTxome),
    fracSecTxome_p0.99 = sum(nReadsWithSecondaryAlignmentsTxome_p0.99) /
      sum(nReadsWithPrimaryAlignmentsTxome_p0.99)) %>%
  dplyr::select(dataset, fracSecGenome, fracSuppGenome,
                fracSecTxome, fracSuppTxome, fracSecTxome_p0.99) %>%
  tidyr::gather(ftype, fracReads, -dataset) %>%
  dplyr::mutate(ref = stringr::str_extract(ftype, "Genome|Txome")) %>%
  dplyr::mutate(ref = replace(ref, ref == "Txome", "Transcriptome")) %>%
  dplyr::mutate(ftype = gsub("Genome|Txome", "", ftype)) %>%
  dplyr::mutate(
    ftype = replace(ftype, ftype == "fracSec", 
                    paste0("Fraction of reads with primary alignment", 
                           " that also have secondary alignment(s)")),
    ftype = replace(ftype, ftype == "fracSupp",
                    paste0("Fraction of reads with primary alignment",
                           " that also have supplementary alignment(s)")),
    ftype = replace(ftype, ftype == "fracSec_p0.99", 
                    paste0("Fraction of reads with primary alignment",
                           " that also have secondary alignment(s)", 
                           ", with -p=0.99"))
  )

colfrac <- c(
  `Fraction of reads with primary alignment that also have secondary alignment(s)` = "#E8601C",
  `Fraction of reads with primary alignment that also have secondary alignment(s), with -p=0.99` = "#7BAFDE",
  `Fraction of reads with primary alignment that also have supplementary alignment(s)` = "#90C987"
)

## Find the fraction of the total number of reads that have secondary alignment
## _and_ where the primary alignment is to a duplicated transcript
dup_tx <- unique(unlist(read.delim(duptx, header = TRUE, as.is = TRUE)))
sec_dup_tx <- 
  dplyr::bind_rows(
    tAlign %>% dplyr::mutate(duplicatedtx = rname %in% dup_tx) %>%
      dplyr::filter(flag %in% c(0, 16)) %>% dplyr::group_by(sample, dataset) %>% 
      dplyr::summarize(fracReads = sum(duplicatedtx & nbrSecondaryAlignments > 0) /
                         length(nbrSecondaryAlignments)) %>%
      dplyr::mutate(ftype = paste0("Fraction of reads with primary alignment ", 
                                   "that also have secondary alignment(s)")) %>%
      dplyr::mutate(ref = "Transcriptome"),
    tAlign_p0.99 %>% dplyr::mutate(duplicatedtx = rname %in% dup_tx) %>%
      dplyr::filter(flag %in% c(0, 16)) %>% dplyr::group_by(sample, dataset) %>% 
      dplyr::summarize(fracReads = sum(duplicatedtx & nbrSecondaryAlignments > 0) /
                         length(nbrSecondaryAlignments)) %>%
      dplyr::mutate(ftype = paste0("Fraction of reads with primary alignment ", 
                                   "that also have secondary alignment(s), ", 
                                   "with -p=0.99")) %>%
      dplyr::mutate(ref = "Transcriptome"),
    nReadsDataFS %>% 
      dplyr::filter(ftype == paste0("Fraction of reads with primary alignment ", 
                                    "that also have supplementary alignment(s)") | 
                      ref == "Genome") %>%
      dplyr::mutate(fracReads = 0)
  )

png(gsub("\\.rds$", "_nsecondary_duplicated.png", outrds),
    width = 10, height = 6, unit = "in", res = 400)
print(ggplot(sec_dup_tx %>% dplyr::filter(fracReads > 0) %>%
               dplyr::mutate(dataset = factor(dataset, levels = 
                                                ds_order[ds_order %in% dataset])), 
             aes(x = sample, y = fracReads, group = ftype)) + 
        geom_bar(stat = "identity", position = "dodge", 
                 aes(fill = ftype, color = ftype)) + 
        scale_y_continuous(expand = c(0, 0, 0.05, 0)) + 
        theme_bw() + facet_grid(ref ~ dataset, scales = "free_x", space = "free_x") + 
        scale_fill_manual(name = "", values = colfrac) + 
        scale_color_manual(name = "", values = colfrac) + 
        theme(legend.position = "bottom",
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
        guides(fill = guide_legend(nrow = 2, byrow = TRUE)) + 
        xlab("") + ylab("") + 
        ggtitle("Fraction of reads with secondary alignment to identical transcript"))
dev.off()


## Find the fraction of the total number of reads that have secondary alignments
## where all the secondary alignments are to the same gene as the primary, or
## all to other genes, or a mix.
f <- function(primary, secondary) {
  if (length(secondary) == 0) {
    return(NA_character_)
  } else if (all(secondary == primary)) {
    return("AllPrimary")
  } else {
    return("Mix")
  }
}
sec_same_diff_gene <- 
  dplyr::bind_rows(
    tAlign %>% dplyr::filter(flag %in% c(0, 16, 256, 272)) %>%
      dplyr::group_by(read, sample, dataset, nbrSecondaryAlignments) %>%
      dplyr::summarize(stype = f(primary = gname[flag %in% c(0, 16)], 
                                 secondary = gname[flag %in% c(256, 272)])) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(sample, dataset) %>% 
      dplyr::mutate(totReads = length(nbrSecondaryAlignments)) %>%
      dplyr::group_by(sample, dataset, stype) %>%
      dplyr::summarize(fracReads = sum(nbrSecondaryAlignments > 0)/totReads[1]) %>%
      dplyr::mutate(ftype = paste0("Fraction of reads with primary alignment ", 
                                   "that also have secondary alignment(s)")) %>%
      dplyr::mutate(ref = "Transcriptome") %>%
      dplyr::filter(stype == "AllPrimary"),
    tAlign_p0.99 %>% dplyr::filter(flag %in% c(0, 16, 256, 272)) %>%
      dplyr::group_by(read, sample, dataset, nbrSecondaryAlignments) %>%
      dplyr::summarize(stype = f(primary = gname[flag %in% c(0, 16)], 
                                 secondary = gname[flag %in% c(256, 272)])) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(sample, dataset) %>% 
      dplyr::mutate(totReads = length(nbrSecondaryAlignments)) %>%
      dplyr::group_by(sample, dataset, stype) %>%
      dplyr::summarize(fracReads = sum(nbrSecondaryAlignments > 0)/totReads[1]) %>%
      dplyr::mutate(ftype = paste0("Fraction of reads with primary alignment ", 
                                   "that also have secondary alignment(s), ", 
                                   "with -p=0.99")) %>%
      dplyr::mutate(ref = "Transcriptome") %>%
      dplyr::filter(stype == "AllPrimary"),
    nReadsDataFS %>% 
      dplyr::filter(ftype == paste0("Fraction of reads with primary alignment ", 
                                    "that also have supplementary alignment(s)") | 
                      ref == "Genome") %>%
      dplyr::mutate(fracReads = 0)
  )

## Same but summarized by dataset
sec_same_diff_gene_ds <- 
  dplyr::bind_rows(
    tAlign %>% dplyr::filter(flag %in% c(0, 16, 256, 272)) %>%
      dplyr::group_by(read, sample, dataset, nbrSecondaryAlignments) %>%
      dplyr::summarize(stype = f(primary = gname[flag %in% c(0, 16)], 
                                 secondary = gname[flag %in% c(256, 272)])) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(dataset) %>% 
      dplyr::mutate(totReads = length(nbrSecondaryAlignments)) %>%
      dplyr::group_by(dataset, stype) %>%
      dplyr::summarize(fracReads = sum(nbrSecondaryAlignments > 0)/totReads[1]) %>%
      dplyr::mutate(ftype = paste0("Fraction of reads with primary alignment ", 
                                   "that also have secondary alignment(s)")) %>%
      dplyr::mutate(ref = "Transcriptome") %>%
      dplyr::filter(stype == "AllPrimary"),
    tAlign_p0.99 %>% dplyr::filter(flag %in% c(0, 16, 256, 272)) %>%
      dplyr::group_by(read, sample, dataset, nbrSecondaryAlignments) %>%
      dplyr::summarize(stype = f(primary = gname[flag %in% c(0, 16)], 
                                 secondary = gname[flag %in% c(256, 272)])) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(dataset) %>% 
      dplyr::mutate(totReads = length(nbrSecondaryAlignments)) %>%
      dplyr::group_by(dataset, stype) %>%
      dplyr::summarize(fracReads = sum(nbrSecondaryAlignments > 0)/totReads[1]) %>%
      dplyr::mutate(ftype = paste0("Fraction of reads with primary alignment ", 
                                   "that also have secondary alignment(s), ", 
                                   "with -p=0.99")) %>%
      dplyr::mutate(ref = "Transcriptome") %>%
      dplyr::filter(stype == "AllPrimary"),
    nReadsDataFSds %>% 
      dplyr::filter(ftype == paste0("Fraction of reads with primary alignment ", 
                                    "that also have supplementary alignment(s)") | 
                      ref == "Genome") %>%
      dplyr::mutate(fracReads = 0)
  )

write.table(sec_same_diff_gene, file = gsub("\\.rds", "_fraction_reads_with_secondary_to_same_gene.txt", outrds), 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(sec_same_diff_gene_ds, file = gsub("\\.rds", "_fraction_reads_with_secondary_to_same_gene_byds.txt", outrds), 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

pfracsupp <- 
  ggplot(nReadsDataFS %>%
           dplyr::ungroup() %>%
           dplyr::mutate(sample = removeDatasetFromSample(sample, dataset)) %>%
           dplyr::mutate(dataset = factor(dataset, levels = ds_order[ds_order %in% dataset])),
         aes(x = sample, y = fracReads, group = ftype)) + 
  geom_bar(stat = "identity", position = "dodge", 
           aes(fill = ftype, color = ftype)) + 
  scale_y_continuous(expand = c(0, 0, 0.05, 0)) + 
  theme_bw() + facet_grid(ref ~ dataset, scales = "free_x", space = "free_x") + 
  scale_fill_manual(name = "", values = colfrac) + 
  scale_color_manual(name = "", values = colfrac) + 
  theme(legend.position = "bottom",
        strip.text = element_text(size = 8)) + 
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) + 
  xlab("") + ylab("Fraction of aligned reads") + 
  geom_bar(data = sec_same_diff_gene %>% dplyr::ungroup() %>%
             dplyr::mutate(sample = removeDatasetFromSample(sample, dataset)) %>%
             dplyr::mutate(dataset = factor(dataset, levels = 
                                              ds_order[ds_order %in% dataset])), 
           color = "white", fill = "white", 
           stat = "identity", position = "dodge") + 
  geom_bar(data = sec_same_diff_gene %>% dplyr::ungroup() %>%
             dplyr::mutate(sample = removeDatasetFromSample(sample, dataset)) %>%
             dplyr::mutate(dataset = factor(dataset, levels = 
                                              ds_order[ds_order %in% dataset])), 
           aes(fill = ftype, color = ftype), 
           stat = "identity", position = "dodge", alpha = 0.3)

pfracsuppds <- 
  ggplot(nReadsDataFSds %>%
           dplyr::ungroup() %>% 
           dplyr::mutate(dataset = factor(dataset, levels = 
                                            ds_order[ds_order %in% dataset])),
         aes(x = dataset, y = fracReads, group = ftype)) + 
  geom_bar(stat = "identity", position = "dodge", 
           aes(fill = ftype, color = ftype)) + 
  scale_y_continuous(expand = c(0, 0, 0.05, 0)) + 
  theme_bw() + facet_grid(ref ~ ., scales = "free_x", space = "free_x") + 
  scale_fill_manual(name = "", values = colfrac) + 
  scale_color_manual(name = "", values = colfrac) + 
  theme(legend.position = "bottom",
        strip.text = element_text(size = 8),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  guides(fill = guide_legend(nrow = 3, byrow = TRUE)) + 
  xlab("") + ylab("Fraction of aligned reads") + 
  geom_bar(data = sec_same_diff_gene_ds %>% dplyr::ungroup() %>%
             dplyr::mutate(dataset = factor(dataset, levels = 
                                              ds_order[ds_order %in% dataset])), 
           color = "white", fill = "white", 
           stat = "identity", position = "dodge") + 
  geom_bar(data = sec_same_diff_gene_ds %>% dplyr::ungroup() %>%
             dplyr::mutate(dataset = factor(dataset, levels = 
                                              ds_order[ds_order %in% dataset])), 
           aes(fill = ftype, color = ftype), 
           stat = "identity", position = "dodge", alpha = 0.3)

## Merge with primary/supplementary pairs
dfprimsupp <- do.call(rbind, lapply(datasets, function(ds) {
  readRDS(paste0(ds, "/output/", ds, 
                 "_primary_supplementary_alignments_distances.rds")) %>% 
    dplyr::mutate(sample = remap[sample]) %>% 
    dplyr::mutate(dataset = remapds[ds])
})) %>%
  dplyr::left_join(sample_annotation %>% dplyr::select(sample_remap, condition),
                   by = c("sample" = "sample_remap")) %>%
  dplyr::filter(condition %in% conditions) %>%
  dplyr::select(-condition) %>%
  dplyr::filter(strands == "same chromosome, different strand" & 
                  distn == 0) %>%
  dplyr::left_join(gAlign %>% dplyr::filter(flag %in% c(0, 16)) %>%
                     dplyr::select(read, sample, dataset, NM, nbrM, nbrS, 
                                   nbrD, nbrI, alignedLength, 
                                   nbrSecondaryAlignments, 
                                   nbrSupplementaryAlignments, fracM, fracS),
                   by = c("sample", "dataset", "read")) %>%
  dplyr::filter(dataset %in% c("ONT-NSK007-HAP", "ONT-DCS108-HAP")) %>%
  dplyr::mutate(dataset = factor(dataset, levels = ds_order[ds_order %in% dataset]))
gps <- ggplot(dfprimsupp, aes(x = dataset, y = ovlap/(nbrM + nbrD),
                              fill = dataset)) + 
  geom_violin() + ylab("Overlap/primary alignment length") + 
  xlab("") + theme_bw() + 
  scale_fill_manual(values = ds_colors) + 
  theme(legend.position = "none")

## ========================================================================== ##
## Read length vs nbr of supplementary alignments
## ========================================================================== ##
pnsupprlbox <- 
  ggplot(gAlign %>% dplyr::filter(flag %in% c(0, 16)) %>%
           dplyr::mutate(nbrSupplementaryAlignments = replace(
             nbrSupplementaryAlignments, nbrSupplementaryAlignments > 10, ">10"
           )) %>%
           dplyr::mutate(nbrSupplementaryAlignments = factor(
             nbrSupplementaryAlignments, levels = c(0:10, ">10")
           )) %>%
           dplyr::mutate(dataset = factor(dataset, levels = 
                                            ds_order[ds_order %in% dataset])) %>%
           dplyr::arrange(dataset, sample) %>%
           dplyr::mutate(sample = factor(sample, levels = unique(sample))),
         aes(x = nbrSupplementaryAlignments, y = readLength, 
             fill = dataset)) + 
  geom_boxplot() + theme_bw() + scale_y_log10() + 
  facet_wrap(~ sample) + 
  scale_fill_manual(values = ds_colors, name = "") + 
  theme(legend.position = "none") + 
  xlab("Number of supplementary alignments") + ylab("Read length")

png(gsub("\\.rds$", "_nsupp_vs_readlength_box.png", outrds),
    width = 12, height = 8.5, unit = "in", res = 400)
print(pnsupprlbox)
dev.off()

png(gsub("\\.rds$", "_readlength_vs_nsupp_density.png", outrds),
    width = 12, height = 8.5, unit = "in", res = 400)
prlvsnsupp <- 
  ggplot(gAlign %>% dplyr::filter(flag %in% c(0, 16)) %>%
           dplyr::mutate(nbrSupplementaryAlignments = replace(
             nbrSupplementaryAlignments, nbrSupplementaryAlignments > 0, ">0"
           )) %>%
           dplyr::mutate(nbrSupplementaryAlignments = factor(
             nbrSupplementaryAlignments, levels = c(0, ">0")
           )) %>%
           dplyr::mutate(dataset = factor(dataset, levels = 
                                            ds_order[ds_order %in% dataset])) %>%
           dplyr::arrange(dataset, sample) %>%
           dplyr::mutate(sample = factor(sample, levels = unique(sample))),
         aes(x = readLength, fill = nbrSupplementaryAlignments)) + 
  geom_density(alpha = 0.5) + theme_bw() + scale_x_log10() + 
  scale_fill_manual(values = c("red", "blue"), 
                    name = "Number of supplementary alignments") + 
  theme(legend.position = "bottom") + 
  xlab("Read length") + ylab("Density")
print(prlvsnsupp + facet_wrap(~ sample))
dev.off()

## ========================================================================== ##
## Aligned length vs number of covered junctions (primary alignments)
## ========================================================================== ##
pnjuncrl <- 
  ggplot(gAlign %>% dplyr::filter(flag %in% c(0, 16)) %>%
           dplyr::mutate(dataset = factor(dataset, levels = 
                                            ds_order[ds_order %in% dataset])) %>%
           dplyr::arrange(dataset, sample) %>%
           dplyr::mutate(sample = factor(sample, levels = unique(sample))),
         aes(x = alignedLength, y = nbrJunctions)) + 
  geom_hex(bins = 100, aes(fill = stat(density))) + theme_bw() + scale_x_log10() + 
  scale_fill_gradient(name = "", low = "bisque2", high = "darkblue") + 
  facet_wrap(~ sample) + 
  xlab("Aligned length") + ylab("Number of spanned junctions")

png(gsub("\\.rds$", "_njunc_vs_alignedlength.png", outrds),
    width = 12, height = 8.5, unit = "in", res = 400)
print(pnjuncrl)
dev.off()

## ========================================================================== ##
## Number of covered junctions (primary alignments)
## ========================================================================== ##
djunc <- gAlign %>% dplyr::filter(flag %in% c(0, 16)) %>%
  dplyr::mutate(nbrJunctions = replace(nbrJunctions, nbrJunctions > 20, ">20")) %>%
  dplyr::mutate(nbrJunctions = factor(nbrJunctions, 
                                      levels = c(as.character(0:20), ">20"))) %>%
  dplyr::group_by(sample, dataset, nbrJunctions) %>%
  dplyr::tally() %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(n = n/sum(n)) %>%
  dplyr::ungroup()

png(gsub("\\.rds$", "_njunc_distribution_bysample.png", outrds),
    width = 12, height = 8.5, unit = "in", res = 400)
print(ggplot(djunc %>%
               dplyr::mutate(dataset = factor(dataset, levels = 
                                                ds_order[ds_order %in% dataset])) %>%
               dplyr::arrange(dataset, sample) %>%
               dplyr::mutate(sample = factor(sample, levels = unique(sample))), 
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

## Combine with Illumina
ilmnnjunc <- readRDS(ilmnnjunc)

png(gsub("\\.rds$", "_njunc_distribution.png", outrds),
    width = 8, height = 3, unit = "in", res = 400)
p1 <- ggplot(djunc, aes(x = nbrJunctions, y = n), 
             fill = "#B3B3B3", color = "#B3B3B3") +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        legend.position = "none") + 
  scale_y_continuous(expand = c(0, 0, 0.05, 0)) + 
  xlab("Number of spanned junctions") + ylab("Fraction of reads") + 
  ggtitle("Nanopore")
p2 <- ggplot(ilmnnjunc, aes(x = nbrJunctions, y = n), 
             fill = "#B3B3B3", color = "#B3B3B3") +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        legend.position = "none") +
  scale_y_continuous(expand = c(0, 0, 0.05, 0)) + 
  xlab("Number of spanned junctions") + ylab("Fraction of reads") + 
  ggtitle("Illumina")
cowplot::plot_grid(p1, p2, nrow = 1, rel_widths = c(1, 1),
                   labels = c("A", "B"), align = "h")
dev.off()

## ========================================================================== ##
## Read length vs aligned length (primary alignments)
## ========================================================================== ##
png(gsub("\\.rds$", "_alignedlength_vs_readlength_genome.png", outrds),
    width = 12, height = 8, unit = "in", res = 400)
print(ggplot(gAlign %>% dplyr::filter(flag %in% c(0, 16)) %>%
               dplyr::mutate(dataset = factor(dataset, levels = 
                                                ds_order[ds_order %in% dataset])) %>%
               dplyr::arrange(dataset, sample) %>%
               dplyr::mutate(sample = factor(sample, levels = unique(sample))),
             aes(x = readLength, y = alignedLength)) + 
        geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
        facet_wrap(~ sample) + geom_hex(bins = 100, aes(fill = stat(density))) + 
        theme_bw() + scale_x_log10() + scale_y_log10() + 
        scale_fill_gradient(name = "", low = "bisque2", high = "darkblue") + 
        xlab("Read length") + ylab("Aligned length (primary alignment)") + 
        ggtitle("Genome"))
dev.off()

png(gsub("\\.rds$", "_alignedlength_vs_readlength_genome_byds_sub.png", outrds),
    width = 6, height = 3, unit = "in", res = 400)
alvsrlgmdssub <- 
  ggplot(gAlign %>% dplyr::filter(flag %in% c(0, 16)) %>%
           dplyr::filter(dataset %in% c("ONT-RNA001-HAP", "ONT-DCS108-HAP")) %>%
           dplyr::mutate(dataset = factor(dataset, levels = 
                                            ds_order[ds_order %in% dataset])),
         aes(x = readLength, y = alignedLength)) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
  facet_wrap(~ dataset) + geom_hex(bins = 100, aes(fill = stat(density))) + 
  theme_bw() + scale_x_log10() + scale_y_log10() + 
  scale_fill_gradient(name = "", low = "bisque2", high = "darkblue") + 
  xlab("Read length") + ylab("Aligned length (primary alignment)") + 
  ggtitle("Genome")
print(alvsrlgmdssub)
dev.off()

png(gsub("\\.rds$", "_alignedlength_vs_readlength_genome_byds_sub_sqrt.png", outrds),
    width = 6, height = 3, unit = "in", res = 400)
alvsrlgmdssubsqrt <- 
  ggplot(gAlign %>% dplyr::filter(flag %in% c(0, 16)) %>%
           dplyr::filter(dataset %in% c("ONT-RNA001-HAP", "ONT-DCS108-HAP")) %>%
           dplyr::mutate(dataset = factor(dataset, levels = 
                                            ds_order[ds_order %in% dataset])),
         aes(x = readLength, y = alignedLength)) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
  facet_wrap(~ dataset) + geom_hex(bins = 100, aes(fill = stat(density))) + 
  theme_bw() + scale_x_sqrt() + scale_y_sqrt() + 
  scale_fill_gradient(name = "", low = "bisque2", high = "darkblue") + 
  xlab("Read length") + ylab("Aligned length (primary alignment)") + 
  ggtitle("Genome")
print(alvsrlgmdssubsqrt)
dev.off()

png(gsub("\\.rds$", "_alignedlength_vs_readlength_genome_reads_not_aligning_to_txome.png", outrds),
    width = 12, height = 8, unit = "in", res = 400)
print(ggplot(gnottAlign %>% dplyr::filter(flag %in% c(0, 16)) %>%
               dplyr::mutate(dataset = factor(dataset, levels = 
                                                ds_order[ds_order %in% dataset])) %>%
               dplyr::arrange(dataset, sample) %>%
               dplyr::mutate(sample = factor(sample, levels = unique(sample))),
             aes(x = readLength, y = alignedLength)) + 
        geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
        facet_wrap(~ sample) + geom_hex(bins = 100, aes(fill = stat(density))) + 
        theme_bw() + scale_x_log10() + scale_y_log10() + 
        scale_fill_gradient(name = "", low = "bisque2", high = "darkblue") + 
        xlab("Read length") + ylab("Aligned length (primary alignment)") + 
        ggtitle("Genome (reads not aligning to the transcriptome)"))
dev.off()

png(gsub("\\.rds$", "_alignedlength_vs_readlength_genome_reads_without_suppl_alignments.png", outrds),
    width = 12, height = 8, unit = "in", res = 400)
alvsrlwosupp <- 
  ggplot(gAlign %>% dplyr::filter(flag %in% c(0, 16)) %>%
           dplyr::filter(nbrSupplementaryAlignments == 0) %>%
           dplyr::mutate(dataset = factor(dataset, levels = 
                                            ds_order[ds_order %in% dataset])) %>%
           dplyr::arrange(dataset, sample) %>%
           dplyr::mutate(sample = factor(sample, levels = unique(sample))),
         aes(x = readLength, y = alignedLength)) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
  facet_wrap(~ sample) + geom_hex(bins = 100, aes(fill = stat(density))) + 
  theme_bw() + scale_x_log10() + scale_y_log10() + 
  theme(strip.text = element_text(size = 6),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  scale_fill_gradient(name = "", low = "bisque2", high = "darkblue") + 
  xlab("Read length") + ylab("Aligned length (primary alignment)") + 
  ggtitle("Genome, reads without supplementary alignments")
print(alvsrlwosupp)
dev.off()

png(gsub("\\.rds$", "_alignedlength_vs_readlength_genome_reads_with_suppl_alignments.png", outrds),
    width = 12, height = 8, unit = "in", res = 400)
alvsrlwsupp <- 
  ggplot(gAlign %>% dplyr::filter(flag %in% c(0, 16)) %>%
           dplyr::filter(nbrSupplementaryAlignments > 0) %>%
           dplyr::mutate(dataset = factor(dataset, levels = 
                                            ds_order[ds_order %in% dataset])) %>%
           dplyr::arrange(dataset, sample) %>%
           dplyr::mutate(sample = factor(sample, levels = unique(sample))),
         aes(x = readLength, y = alignedLength)) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
  facet_wrap(~ sample) + geom_hex(bins = 100, aes(fill = stat(density))) + 
  theme_bw() + scale_x_log10() + scale_y_log10() + 
  theme(strip.text = element_text(size = 6),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  scale_fill_gradient(name = "", low = "bisque2", high = "darkblue") + 
  xlab("Read length") + ylab("Aligned length (primary alignment)") + 
  ggtitle("Genome, reads with supplementary alignments")
print(alvsrlwsupp)
dev.off()

png(gsub("\\.rds$", "_alignedlength_vs_readlength_txome.png", outrds),
    width = 12, height = 8, unit = "in", res = 400)
print(ggplot(tAlign_p0.99 %>% dplyr::filter(flag %in% c(0, 16)) %>%
               dplyr::mutate(dataset = factor(dataset, levels = 
                                                ds_order[ds_order %in% dataset])) %>%
               dplyr::arrange(dataset, sample) %>%
               dplyr::mutate(sample = factor(sample, levels = unique(sample))),
             aes(x = readLength, y = alignedLength)) + 
        geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
        facet_wrap(~ sample) + geom_hex(bins = 100, aes(fill = stat(density))) + 
        theme_bw() + scale_x_log10() + scale_y_log10() + 
        scale_fill_gradient(name = "", low = "bisque2", high = "darkblue") + 
        xlab("Read length") + ylab("Aligned length (primary alignment)") + 
        ggtitle("Transcriptome")
)
dev.off()

## Combine with read length vs nbr suppl, and Illumina coverage degree
txlvsbc <- readRDS(txlvsbc)
png(gsub("\\.rds$", "_alignedlength_vs_readlength_vs_nbrsuppl_vs_covdegree.png", outrds),
    width = 12, height = 12, unit = "in", res = 400)
print(
  cowplot::plot_grid(
    cowplot::plot_grid(
      alvsrlwosupp,
      alvsrlwsupp,
      nrow = 1, rel_widths = c(1, 1), labels = c("A", "B")
    ),
    cowplot::plot_grid(
      prlvsnsupp + facet_wrap(~ dataset),
      txlvsbc$txlength_vs_basecoverage_minTPM1,
      nrow = 1, rel_widths = c(1, 1), labels = c("C", "D"),
      align = "h", axis = "tb"
    ),
    ncol = 1, rel_heights = c(1, 1)
  )
)
dev.off()

## Combine with nbr aligned, frac sec/supp, prim/supp overlap
primsupp <- readRDS(primsupp)
png(gsub("\\.rds$", "_naligned_reads_primary_supplementary.png", outrds),
    width = 15, height = 12, unit = "in", res = 400)
print(
  cowplot::plot_grid(pnreads, pfracsupp, primsupp$primary_supplementary_fraction,
                     ncol = 1, rel_heights = c(1, 1, 1), labels = c("A", "B", "C"),
                     align = "v", axis = "lr")
)
dev.off()

png(gsub("\\.rds$", "_naligned_reads_primary_supplementary_byds.png", outrds),
    width = 17, height = 11, unit = "in", res = 400)
print(
  cowplot::plot_grid(
    cowplot::plot_grid(
      pnreadsds + 
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1)), 
      pfracsuppds + 
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1)), 
      primsupp$primary_supplementary_fraction_byds + 
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1)),
      nrow = 1, rel_widths = c(1, 1, 1), labels = c("A", "B", "C"),
      align = "h", axis = "tb"),
    cowplot::plot_grid(
      gps,
      alvsrlgmdssub,
      nrow = 1, rel_widths = c(1, 2), labels = c("D", "E"), align = "h", axis = "tb"
    ),
    ncol = 1, rel_heights = c(1.4, 1), labels = c("", ""), align = "v", axis = "l"
  )
)
dev.off()

png(gsub("\\.rds$", "_naligned_reads_primary_supplementary_byds_sqrt.png", outrds),
    width = 17, height = 11, unit = "in", res = 400)
print(
  cowplot::plot_grid(
    cowplot::plot_grid(
      pnreadsds + 
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1)), 
      pfracsuppds + 
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1)), 
      primsupp$primary_supplementary_fraction_byds + 
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1)),
      nrow = 1, rel_widths = c(1, 1, 1), labels = c("A", "B", "C"),
      align = "h", axis = "tb"),
    cowplot::plot_grid(
      gps,
      alvsrlgmdssubsqrt,
      nrow = 1, rel_widths = c(1, 2), labels = c("D", "E"), align = "h", axis = "tb"
    ),
    ncol = 1, rel_heights = c(1.4, 1), labels = c("", ""), align = "v", axis = "l"
  )
)
dev.off()

## ========================================================================== ##
## Read length vs transcript length (primary alignments)
## ========================================================================== ##
png(gsub("\\.rds$", "_txlength_vs_readlength_txome.png", outrds),
    width = 12, height = 8, unit = "in", res = 400)
print(ggplot(tAlign_p0.99 %>% dplyr::filter(flag %in% c(0, 16)) %>%
               dplyr::mutate(dataset = factor(dataset, levels = 
                                                ds_order[ds_order %in% dataset])) %>%
               dplyr::arrange(dataset, sample) %>%
               dplyr::mutate(sample = factor(sample, levels = unique(sample))),
             aes(x = txLength, y = readLength)) + 
        geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
        facet_wrap(~ sample) + 
        geom_hex(bins = 100, aes(fill = stat(density))) + 
        theme_bw() + scale_x_log10() + scale_y_log10() + 
        scale_fill_gradient(name = "", low = "bisque2", high = "darkblue") + 
        xlab("Transcript length") + ylab("Read length") + ggtitle("Transcriptome")
)
dev.off()

## ========================================================================== ##
## Nbr soft-clipped bases vs transcript length (primary alignments)
## ========================================================================== ##
png(gsub("\\.rds$", "_txlength_vs_S_txome.png", outrds),
    width = 12, height = 8, unit = "in", res = 400)
print(ggplot(tAlign_p0.99 %>% dplyr::filter(flag %in% c(0, 16)) %>%
               dplyr::mutate(dataset = factor(dataset, levels = 
                                                ds_order[ds_order %in% dataset])) %>%
               dplyr::arrange(dataset, sample) %>%
               dplyr::mutate(sample = factor(sample, levels = unique(sample))),
             aes(x = txLength, y = fracS)) + 
        facet_wrap(~ sample) + geom_hex(bins = 100, aes(fill = stat(density))) + 
        theme_bw() + scale_x_log10() + 
        scale_fill_gradient(name = "", low = "bisque2", high = "darkblue") + 
        xlab("Transcript length, primary alignment") + 
        ylab("Fraction soft-clipped bases") + 
        ggtitle("Transcriptome")
)
dev.off()

png(gsub("\\.rds$", "_txlength_vs_S_txome_log10.png", outrds),
    width = 12, height = 8, unit = "in", res = 400)
print(ggplot(tAlign_p0.99 %>% dplyr::filter(flag %in% c(0, 16)) %>%
               dplyr::mutate(dataset = factor(dataset, levels = 
                                                ds_order[ds_order %in% dataset])) %>%
               dplyr::arrange(dataset, sample) %>%
               dplyr::mutate(sample = factor(sample, levels = unique(sample))),
             aes(x = txLength, y = fracS)) + 
        facet_wrap(~ sample) + 
        geom_hex(bins = 100, aes(fill = stat(density))) + 
        theme_bw() + scale_x_log10() + scale_y_log10() + 
        scale_fill_gradient(name = "", low = "bisque2", high = "darkblue") + 
        xlab("Transcript length, primary alignment") + 
        ylab("Fraction soft-clipped bases") + 
        ggtitle("Transcriptome")
)
dev.off()

## ========================================================================== ##
## Average base quality, read length, frac M, frac S distributions
## ========================================================================== ##
df <- dplyr::bind_rows(
  allReads %>% dplyr::select(sample, dataset, rtype, readLength, aveBaseQuality),
  gAlign %>% dplyr::filter(flag %in% c(0, 16)) %>% 
    dplyr::left_join(allReads %>% dplyr::select(sample, read, aveBaseQuality),
                     by = c("sample", "read")) %>%
    dplyr::select(sample, dataset, rtype, readLength, aveBaseQuality,
                  fracM, fracS, fracI, accuracy),
  tAlign_p0.99 %>% dplyr::filter(flag %in% c(0, 16)) %>% 
    dplyr::left_join(allReads %>% dplyr::select(sample, read, aveBaseQuality),
                     by = c("sample", "read")) %>%
    dplyr::select(sample, dataset, rtype, readLength, aveBaseQuality,
                  fracM, fracS, fracI, accuracy),
  allReads %>% dplyr::filter(!(read %in% c(gAlign$read, tAlign_p0.99$read))) %>%
    dplyr::select(sample, dataset, rtype, readLength, aveBaseQuality) %>%
    dplyr::mutate(rtype = "Unaligned reads"),
  gnottAlign %>% dplyr::filter(flag %in% c(0, 16)) %>%
    dplyr::left_join(allReads %>% dplyr::select(sample, read, aveBaseQuality),
                     by = c("sample", "read")) %>%
    dplyr::select(sample, dataset, rtype, readLength, aveBaseQuality,
                  fracM, fracS, fracI, accuracy)
) %>%
  dplyr::mutate(dataset = factor(dataset, levels = ds_order[ds_order %in% dataset]))

colread <- c(
  `All reads` = "#777777", 
  `Reads aligning to the genome` = "#E8601C", 
  `Reads aligning to the transcriptome` = "#7BAFDE",
  `Reads aligning to the transcriptome (-p 0.99)` = "#90C987",
  `Unaligned reads` = "lightblue",
  `Reads aligning to the genome but not the transcriptome` = "#dd99ff"
)

pbqd <- ggplot(df, aes(x = aveBaseQuality, group = interaction(sample, rtype),
                       color = rtype)) + 
  geom_line(stat = "density", size = 1) + 
  facet_wrap(~ dataset, nrow = 1) + theme_bw() + 
  scale_color_manual(values = colread, name = "") + 
  xlab("Average base quality") + ylab("Density") + 
  theme(legend.position = "bottom")

prld <- ggplot(df, aes(x = readLength, group = interaction(sample, rtype), 
                       color = rtype)) + 
  geom_line(stat = "density", size = 1) + scale_x_log10() + 
  facet_wrap(~ dataset, nrow = 1) + theme_bw() + 
  scale_color_manual(values = colread, name = "") + 
  xlab("Read length") + ylab("Density") + 
  theme(legend.position = "bottom")

prldsqrt <- ggplot(df, aes(x = readLength, group = interaction(sample, rtype), 
                           color = rtype)) + 
  geom_line(stat = "density", size = 1) + scale_x_sqrt() + 
  facet_wrap(~ dataset, nrow = 1) + theme_bw() + 
  scale_color_manual(values = colread, name = "") + 
  xlab("Read length") + ylab("Density") + 
  theme(legend.position = "bottom")

pmd <- ggplot(df, aes(x = fracM, group = interaction(sample, rtype), color = rtype)) + 
  geom_line(stat = "density", size = 1) + 
  facet_wrap(~ dataset, nrow = 1) + theme_bw() + 
  scale_color_manual(values = colread, name = "") + 
  xlab("Fraction Ms in primary alignments") + ylab("Density") + 
  theme(legend.position = "bottom")

psd <- ggplot(df, aes(x = fracS, group = interaction(sample, rtype), color = rtype)) + 
  geom_line(stat = "density", size = 1) + 
  facet_wrap(~ dataset, nrow = 1) + theme_bw() + 
  scale_color_manual(values = colread, name = "") + 
  xlab("Fraction soft-clipped bases in primary alignments") + ylab("Density") + 
  theme(legend.position = "bottom")

pid <- ggplot(df, aes(x = fracI, group = interaction(sample, rtype), color = rtype)) + 
  geom_line(stat = "density", size = 1) + 
  facet_wrap(~ dataset, nrow = 1) + theme_bw() + 
  scale_color_manual(values = colread, name = "") + 
  xlab("Fraction insertions in primary alignments") + ylab("Density") + 
  theme(legend.position = "bottom")

pacc <- ggplot(df, aes(x = accuracy, group = interaction(sample, rtype), color = rtype)) + 
  geom_line(stat = "density", size = 1) + 
  facet_wrap(~ dataset, nrow = 1) + theme_bw() + 
  scale_color_manual(values = colread, name = "") + 
  xlab("Base accuracy in primary alignments") + ylab("Density") + 
  theme(legend.position = "bottom")

png(gsub("\\.rds$", "_read_length_and_MSI_rate.png", outrds), 
    width = 12, height = 15, unit = "in", res = 400)
print(cowplot::plot_grid(pbqd, prld, pmd, psd, pacc, 
                         ncol = 1, labels = c("A", "B", "C", "D", "E")))
dev.off()

png(gsub("\\.rds$", "_read_length_and_MSI_rate_sqrt.png", outrds), 
    width = 12, height = 15, unit = "in", res = 400)
print(cowplot::plot_grid(pbqd, prldsqrt, pmd, psd, pacc, 
                         ncol = 1, labels = c("A", "B", "C", "D", "E")))
dev.off()

png(gsub("\\.rds$", "_avebasequality_vs_accuracy_primary_genome_alignments.png", outrds), 
    width = 12, height = 12, unit = "in", res = 400)
print(ggplot(df %>% dplyr::filter(rtype == "Reads aligning to the genome") %>%
               dplyr::arrange(dataset, sample) %>%
               dplyr::mutate(sample = factor(sample, levels = unique(sample))), 
             aes(x = aveBaseQuality, y = accuracy)) + 
        geom_hex(bins = 100, aes(fill = stat(density))) + 
        scale_fill_gradient(name = "", low = "bisque2", high = "darkblue") + 
        facet_wrap(~ sample, scales = "free") + theme_bw() + 
        xlab("Average base quality") + ylab("Base accuracy") + 
        ggtitle("Primary alignments to the genome"))
dev.off()

png(gsub("\\.rds$", "_readlength_vs_accuracy_primary_genome_alignments.png", outrds), 
    width = 12, height = 12, unit = "in", res = 400)
print(ggplot(df %>% dplyr::filter(rtype == "Reads aligning to the genome") %>%
               dplyr::arrange(dataset, sample) %>%
               dplyr::mutate(sample = factor(sample, levels = unique(sample))), 
             aes(x = readLength, y = accuracy)) + 
        geom_hex(bins = 100, aes(fill = stat(density))) + scale_x_log10() + 
        scale_fill_gradient(name = "", low = "bisque2", high = "darkblue") + 
        facet_wrap(~ sample) + theme_bw() + 
        xlab("Read length") + ylab("Base accuracy") + 
        ggtitle("Primary alignments to the genome"))
dev.off()

png(gsub("\\.rds$", "_readlength_vs_accuracy_primary_genome_alignments_sqrt.png", outrds), 
    width = 12, height = 12, unit = "in", res = 400)
print(ggplot(df %>% dplyr::filter(rtype == "Reads aligning to the genome") %>%
               dplyr::arrange(dataset, sample) %>%
               dplyr::mutate(sample = factor(sample, levels = unique(sample))), 
             aes(x = readLength, y = accuracy)) + 
        geom_hex(bins = 100, aes(fill = stat(density))) + scale_x_sqrt() + 
        scale_fill_gradient(name = "", low = "bisque2", high = "darkblue") + 
        facet_wrap(~ sample) + theme_bw() + 
        xlab("Read length") + ylab("Base accuracy") + 
        ggtitle("Primary alignments to the genome"))
dev.off()

## ========================================================================== ##
## Number of supplementary/secondary reads
## ========================================================================== ##
png(gsub("\\.rds$", "_nbr_secondary_genome_alignments.png", outrds),
    width = 10, height = 7, unit = "in", res = 400)
print(ggplot(
  gAlign %>% dplyr::filter(flag %in% c(0, 16)) %>%
    dplyr::mutate(dataset = factor(dataset, levels = ds_order[ds_order %in% dataset])),
  aes(x = nbrSecondaryAlignments, y = stat(density), fill = dataset)) +
    geom_histogram(bins = 11) + 
    theme_bw() + theme(legend.position = "none", 
                       strip.text = element_text(size = 9)) + 
    scale_fill_manual(values = ds_colors, name = "") + 
    xlab("Number of secondary genomic alignments") + 
    facet_wrap(~ sample))
dev.off()

png(gsub("\\.rds$", "_nbr_secondary_transcriptome_alignments.png", outrds),
    width = 10, height = 7, unit = "in", res = 400)
print(ggplot(
  tAlign %>% dplyr::filter(flag %in% c(0, 16)) %>%
    dplyr::mutate(dataset = factor(dataset, levels = ds_order[ds_order %in% dataset])),
  aes(x = nbrSecondaryAlignments, y = stat(density), fill = dataset)) +
    geom_histogram(bins = 25) + 
    theme_bw() + theme(legend.position = "none", 
                       strip.text = element_text(size = 9)) + 
    scale_fill_manual(values = ds_colors, name = "") + 
    xlab("Number of secondary transcriptomic alignments") + 
    facet_wrap(~ sample))
dev.off()

png(gsub("\\.rds$", "_nbr_secondary_transcriptome_alignments_p0.99.png", outrds),
    width = 10, height = 7, unit = "in", res = 400)
print(ggplot(
  tAlign_p0.99 %>% dplyr::filter(flag %in% c(0, 16)) %>%
    dplyr::mutate(dataset = factor(dataset, levels = ds_order[ds_order %in% dataset])),
  aes(x = nbrSecondaryAlignments, y = stat(density), fill = dataset)) +
    geom_histogram(bins = 25) + 
    theme_bw() + theme(legend.position = "none", 
                       strip.text = element_text(size = 9)) + 
    scale_fill_manual(values = ds_colors, name = "") + 
    xlab("Number of secondary transcriptomic alignments (-p=0.99)") + 
    facet_wrap(~ sample))
dev.off()

png(gsub("\\.rds$", "_nbr_secondary_transcriptome_alignments_p0.99_line.png", outrds),
    width = 9, height = 7, unit = "in", res = 400)
print(ggplot(
  tAlign_p0.99 %>% dplyr::filter(flag %in% c(0, 16)) %>%
    dplyr::group_by(sample, dataset, nbrSecondaryAlignments) %>%
    dplyr::tally() %>%
    dplyr::ungroup() %>%
    dplyr::group_by(sample, dataset) %>%
    dplyr::mutate(fraction = n/sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(dataset = factor(dataset, levels = ds_order[ds_order %in% dataset])),
  aes(x = nbrSecondaryAlignments, y = fraction, color = dataset, group = sample)) +
    geom_line() + 
    theme_bw() + 
    scale_color_manual(values = ds_colors, name = "") + 
    xlab("Number of secondary transcriptomic alignments (-p=0.99)") + 
    ylab("Fraction of aligned reads"))
dev.off()

png(gsub("\\.rds$", "_nbr_secondary_transcriptome_alignments_p0.99_zoom.png", outrds),
    width = 10, height = 7, unit = "in", res = 400)
print(ggplot(
  tAlign_p0.99 %>% dplyr::filter(flag %in% c(0, 16)) %>%
    dplyr::mutate(dataset = factor(dataset, levels = ds_order[ds_order %in% dataset])),
  aes(x = nbrSecondaryAlignments, y = stat(density), fill = dataset)) +
    geom_histogram(bins = 101) + coord_cartesian(xlim = c(0, 20)) +  
    theme_bw() + theme(legend.position = "none", 
                       strip.text = element_text(size = 9)) + 
    scale_fill_manual(values = ds_colors, name = "") + 
    xlab("Number of secondary transcriptomic alignments (-p=0.99)") + 
    facet_wrap(~ sample))
dev.off()


png(gsub("\\.rds$", "_nbr_secondary_genome_transcriptome_alignments_p0.99_zoom_line.png", outrds),
    width = 10, height = 3, unit = "in", res = 400)
p1 <- ggplot(
  gAlign %>% dplyr::filter(flag %in% c(0, 16)) %>%
    dplyr::group_by(sample, dataset, nbrSecondaryAlignments) %>%
    dplyr::tally() %>%
    dplyr::ungroup() %>%
    dplyr::group_by(sample, dataset) %>%
    dplyr::mutate(fraction = n/sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(dataset = factor(dataset, levels = ds_order[ds_order %in% dataset])),
  aes(x = nbrSecondaryAlignments, y = fraction, color = dataset, group = sample)) +
  geom_line(size = 1.25) + 
  theme_bw() + 
  scale_color_manual(values = ds_colors, name = "") + 
  xlab("Number of secondary genomic alignments") + 
  ylab("Fraction of aligned reads")
p2 <- ggplot(
  tAlign_p0.99 %>% dplyr::filter(flag %in% c(0, 16)) %>%
    dplyr::group_by(sample, dataset, nbrSecondaryAlignments) %>%
    dplyr::tally() %>%
    dplyr::ungroup() %>%
    dplyr::group_by(sample, dataset) %>%
    dplyr::mutate(fraction = n/sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(dataset = factor(dataset, levels = ds_order[ds_order %in% dataset])),
  aes(x = nbrSecondaryAlignments, y = fraction, color = dataset, group = sample)) +
    geom_line(size = 1.25) + coord_cartesian(xlim = c(0, 20)) +  
    theme_bw() + 
    scale_color_manual(values = ds_colors, name = "") + 
    xlab("Number of secondary transcriptomic alignments (-p=0.99)") + 
    ylab("Fraction of aligned reads")
cowplot::plot_grid(p1 + theme(legend.position = "none"), 
                   p2 + theme(legend.position = "none"),
                   cowplot::get_legend(p2),
                   labels = c("A", "B", ""), nrow = 1, 
                   rel_widths = c(1, 1, 0.4))
dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()




args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(ggbeeswarm)
})

print(tx2gene)
print(outrds)

## -------------------------------------------------------------------------- ##
## Read tx2gene and extract some categories
tx2gene <- readRDS(tx2gene)
tx2gene <- tx2gene %>%
  dplyr::mutate(tx_biotype = replace(tx_biotype, tx_biotype == "protein_coding" & 
                                       chromosome == "MT", "Mt_protein_coding")) %>%
  dplyr::mutate(tx_biotype = replace(tx_biotype, grepl("^RPL[0-9]|^RPS[0-9]", 
                                                       symbol), "ribosomal_protein"))

## -------------------------------------------------------------------------- ##
## Read polyA tail length estimates
tailfindr <- read.delim("RNA001/tailfindr/wt_1_RNA001_tailfindr.csv",
                        header = TRUE, as.is = TRUE, sep = ",")
nanopolish <- read.delim(paste0("RNA001/nanopolish/wt_1_RNA001/", 
                                "pipeline-polya-ng/tails/filtered_tails.tsv"), 
                         header = TRUE, as.is = TRUE)

table(nanopolish$polya_length < 1000)

## -------------------------------------------------------------------------- ##
## Add transcript annotation to nanopolish table
nanopolish <- nanopolish %>%
  dplyr::left_join(tx2gene, by = c("contig" = "tx")) %>%
  dplyr::filter(polya_length < 1000)

## -------------------------------------------------------------------------- ##
## Merge nanopolish and tailfindr estimates
tails <- dplyr::inner_join(
  tailfindr %>% dplyr::select(read_id, tail_length) %>%
    dplyr::rename(readname = read_id, tailfindr = tail_length),
  nanopolish %>% dplyr::select(readname, polya_length) %>%
    dplyr::rename(Nanopolish = polya_length),
  by = "readname"
) %>%
  dplyr::filter(!is.na(tailfindr) & !is.na(Nanopolish))

## -------------------------------------------------------------------------- ##
## Summary plot
png(gsub("\\.rds$", "_summary.png", outrds),
    width = 7, height = 9, unit = "in", res = 400)
cowplot::plot_grid(
  cowplot::plot_grid(
    ggplot(tails %>% 
             tidyr::gather(key = "method", value = "polya_length", -readname), 
           aes(x = method, y = polya_length)) + 
      geom_boxplot() + theme_bw() + 
      xlab("") + ylab("polyA tail length estimate"),
    ggplot(tails, aes(x = Nanopolish, y = tailfindr)) + 
      geom_abline(slope = 1, intercept = 0) + 
      geom_point(alpha = 0.25) + theme_bw() + 
      xlab("Nanopolish polyA tail length estimate") + 
      ylab("tailfindr polyA tail length estimate"),
    rel_widths = c(0.5, 1), nrow = 1, labels = c("A", "B"),
    align = "h", axis = "bt"),
  ggplot(nanopolish %>% dplyr::group_by(tx_biotype) %>%
           dplyr::mutate(nreads = length(readname)) %>%
           dplyr::filter(nreads >= 10), 
         aes(x = tx_biotype, y = polya_length)) + 
    geom_boxplot() + theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    xlab("Transcript biotype") + 
    ylab("Nanopolish polyA tail length estimate"),
  rel_heights = c(1, 1.5), ncol = 1, labels = c("", "C")
)
dev.off()

## -------------------------------------------------------------------------- ##
## Calculate median polyA tail length estimate for reads assigned to each
## transcript. Plot these values for genes with many annotated transcripts. 
np <- nanopolish %>% dplyr::group_by(contig, symbol) %>%
  dplyr::summarize(median_polya_length = median(polya_length)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(symbol) %>%
  dplyr::mutate(ntx = length(contig)) %>%
  dplyr::filter(ntx > 1)

png(gsub("\\.rds$", "_polya_length_divergence_within_gene.png", outrds),
    height = 4, width = 7, unit = "in", res = 400)
ggplot(np %>% dplyr::filter(ntx >= 7), 
       aes(x = symbol, y = median_polya_length, color = symbol)) + 
  geom_beeswarm(size = 2.5, alpha = 0.7) + theme_bw() + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  ylab("Median polyA tail length estimate\nfor reads assigned to transcript") + 
  xlab("Gene")
dev.off()
## -------------------------------------------------------------------------- ##

saveRDS(NULL, file = outrds)
date()
sessionInfo()





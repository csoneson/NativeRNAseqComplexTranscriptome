args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
})
source("manuscript_results_figures/Rscripts/remap_sample_names.R")

datasets <- strsplit(datasets, ",")[[1]]
datasets <- c(datasets, "Illumina")
conditions <- strsplit(conditions, ",")[[1]]

print(datasets)
print(conditions)
print(outrds)

salmonfiles <- unlist(lapply(datasets, function(ds) {
  lf <- list.files(file.path(ds, "salmon31"), pattern = "eq_classes.txt",
                   full.names = TRUE, recursive = TRUE)
  names(lf) <- remap[basename(gsub("/aux_info/eq_classes.txt$", "", lf))]
  lf
}))
salmonfiles <- salmonfiles[names(salmonfiles) %in% 
                             sample_annotation$sample_remap[sample_annotation$condition %in% conditions]]
salmonfiles

eqcl <- lapply(salmonfiles, function(f) {
  x <- readLines(f)
  n_tr <- as.numeric(x[1])  ## Total number of transcripts
  n_eq <- as.numeric(x[2])  ## Total number of equivalence classes
  tx_id <- x[3:(n_tr + 2)]  ## Transcript IDs
  quants <- x[(n_tr + 3):length(x)]  ## Characteristics of equivalence classes
  
  ## Split equivalence class characteristics. Each element of the list corresponds
  ## to one equivalence class, and lists its number of transcripts, the
  ## transcripts IDs and the total number of reads
  do.call(dplyr::bind_rows, lapply(quants, function(w) {
    tmp = strsplit(w, "\\\t")[[1]]
    nbr_tx = as.numeric(tmp[1])
    data.frame(nbr_tx = rep(nbr_tx, as.numeric(tmp[length(tmp)])),
               stringsAsFactors = FALSE)
  }))
})

for (nm in names(eqcl)) {
  eqcl[[nm]]$sample <- nm
  eqcl[[nm]]$dataset <- strsplit(nm, "_")[[1]][1]
}

eqcl <- do.call(dplyr::bind_rows, eqcl) %>%
  dplyr::mutate(dataset = factor(dataset, levels = ds_order[ds_order %in% dataset]))

png(gsub("\\.rds$", ".png", outrds), height = 12, width = 16, 
    unit = "in", res = 400)
p1 <- ggplot(eqcl %>% 
               dplyr::mutate(sample = removeDatasetFromSample(sample, dataset)), 
             aes(x = sample, y = nbr_tx, fill = dataset)) + 
  geom_boxplot() + theme_bw() + xlab("") + 
  ylab("Number of transcripts in equivalence class") + 
  theme(legend.position = "none") + 
  facet_grid(~ dataset, scales = "free_x", space = "free_x") + 
  scale_fill_manual(values = ds_colors, name = "")
p2 <- p1 + 
  stat_summary(fun.y = mean, geom = "point", shape = 18, 
               size = 4, color = "black", fill = "black") +
  coord_cartesian(ylim = c(0, 15))
cowplot::plot_grid(p1 + ggtitle("Full range"), 
                   p2 + ggtitle("Zoomed in"), 
                   ncol = 1, labels = c("A", "B"), rel_heights = c(1, 1))
dev.off()

## By dataset
p1ds <- ggplot(eqcl, 
               aes(x = dataset, y = nbr_tx, fill = dataset)) + 
  geom_boxplot() + theme_bw() + xlab("") + 
  ylab("Number of transcripts in equivalence class") + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  scale_fill_manual(values = ds_colors, name = "")
p2ds <- p1ds + 
  stat_summary(fun.y = mean, geom = "point", shape = 18, 
               size = 4, color = "black", fill = "black") +
  coord_cartesian(ylim = c(0, 15))


saveRDS(list(ptxeqfull = p1 + ggtitle("Full range"), 
             ptxeqzoom = p2 + ggtitle("Zoomed in"),
             ptxeqfullds = p1ds + ggtitle("Full range"), 
             ptxeqzoomds = p2ds + ggtitle("Zoomed in")), file = outrds)
date()
sessionInfo()
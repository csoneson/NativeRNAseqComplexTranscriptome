args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(BSgenome.Hsapiens.NCBI.GRCh38, lib.loc = "/home/charlotte/R/x86_64-pc-linux-gnu-library/3.5")
  library(jcc, lib.loc = "/home/charlotte/R/x86_64-pc-linux-gnu-library/3.5")
  library(dplyr)
})

print(gtf)
print(bam)
print(genomeVersion)
print(version)
print(tx2gene)
print(nCores)
print(predcovrds)

## Fit fragment bias model
biasMod <- fitAlpineBiasModel(gtf = gtf, bam = bam, organism = "Homo_sapiens",
                              genome = Hsapiens, genomeVersion = genomeVersion,
                              version = version, minLength = 600, maxLength = 7000,
                              minCount = 500, maxCount = 10000, 
                              subsample = TRUE, nbrSubsample = 100, seed = 1,
                              minSize = NULL, maxSize = NULL, 
                              verbose = TRUE)

tx2gene <- readRDS(tx2gene)

## Predict transcript coverage profiles
predCovProfiles <- predictTxCoverage(biasModel = biasMod$biasModel, 
                                     exonsByTx = biasMod$exonsByTx, 
                                     bam = bam, tx2gene = tx2gene, 
                                     genome = Hsapiens, genes = NULL, 
                                     nCores = nCores, verbose = TRUE)

saveRDS(predCovProfiles, file = predcovrds)
date()
sessionInfo()

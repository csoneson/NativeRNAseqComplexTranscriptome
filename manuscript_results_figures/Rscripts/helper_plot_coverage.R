
#' Plot JCC score summary for a given gene
#'
#' This function generates a multi-panel plot, including (i) the observed
#' coverage and gene model, (ii) the estimated transcript abundances and (iii)
#' the observed and predicted junction coverages within the gene.
#'
#' @param gtf Path to gtf file with genomic features. Preferably in Ensembl
#'   format.
#' @param junctionCovs A \code{data.frame} with junction coverage information,
#'   output from \code{\link{calculateJCCScores}}.
#' @param useGene A character string giving the name of the gene to plot.
#' @param bwFiles Path to one or more bigwig files generated from a genome
#'   alignment BAM file.
#' @param txQuants A \code{data.frame} with estimated transcript abundances,
#'   output from \code{\link{combineCoverages}}.
#' @param txCol The column in the gtf file that contains the transcript ID.
#' @param geneCol The column in the gtf file that contains the gene ID.
#' @param exonCol The column in the gtf file that contains the exon ID.
#'
#' @return A \code{ggplot} object with a multi-panel plot. The size is optimized
#'   for display in a device of size approximately 12in wide and 10in high.
#'
#' @author Charlotte Soneson
#'
#' @export
#'
#' @references
#' Soneson C, Love MI, Patro R, Hussain S, Malhotra D, Robinson MD: A junction
#' coverage compatibility score to quantify the reliability of transcript
#' abundance estimates and annotation catalogs. bioRxiv doi:10.1101/378539
#' (2018)
#'
#' @examples
#' \dontrun{
#' gtf <- system.file("extdata/Homo_sapiens.GRCh38.90.chr22.gtf.gz",
#'                    package = "jcc")
#' bam <- system.file("extdata/reads.chr22.bam", package = "jcc")
#' biasMod <- fitAlpineBiasModel(gtf = gtf, bam = bam,
#'                               organism = "Homo_sapiens",
#'                               genome = Hsapiens, genomeVersion = "GRCh38",
#'                               version = 90, minLength = 230,
#'                               maxLength = 7000, minCount = 10,
#'                               maxCount = 10000, subsample = TRUE,
#'                               nbrSubsample = 30, seed = 1, minSize = NULL,
#'                               maxSize = 220, verbose = TRUE)
#' tx2gene <- readRDS(system.file("extdata/tx2gene.sub.rds", package = "jcc"))
#' predCovProfiles <- predictTxCoverage(biasModel = biasMod$biasModel,
#'                                      exonsByTx = biasMod$exonsByTx,
#'                                      bam = bam, tx2gene = tx2gene,
#'                                      genome = Hsapiens,
#'                                      genes = c("ENSG00000070371",
#'                                                "ENSG00000093010"),
#'                                      nCores = 1, verbose = TRUE)
#' txQuants <- readRDS(system.file("extdata/quant.sub.rds", package = "jcc"))
#' txsc <- scaleTxCoverages(txCoverageProfiles = predCovProfiles,
#'                          txQuants = txQuants, tx2gene = tx2gene,
#'                          strandSpecific = TRUE, methodName = "Salmon",
#'                          verbose = TRUE)
#' jcov <- read.delim(system.file("extdata/sub.SJ.out.tab", package = "jcc"),
#'                    header = FALSE, as.is = TRUE) %>%
#'  setNames(c("seqnames", "start", "end", "strand", "motif", "annot",
#'            "uniqreads", "mmreads", "maxoverhang")) %>%
#'  dplyr::mutate(strand = replace(strand, strand == 1, "+")) %>%
#'  dplyr::mutate(strand = replace(strand, strand == 2, "-")) %>%
#'  dplyr::select(seqnames, start, end, strand, uniqreads, mmreads) %>%
#'  dplyr::mutate(seqnames = as.character(seqnames))
#' combCov <- combineCoverages(junctionCounts = jcov,
#'                             junctionPredCovs = txsc$junctionPredCovs,
#'                             txQuants = txsc$txQuants)
#' jcc <- calculateJCCScores(junctionCovs = combCov$junctionCovs,
#'                           geneQuants = combCov$geneQuants)
#' bwFile <- system.file("extdata/reads.chr22.bw", package = "jcc")
#' plotGeneSummary(gtf = gtf, junctionCovs = jcc$junctionCovs,
#'                 useGene = "ENSG00000070371", bwFile = bwFile,
#'                 txQuants = combCov$txQuants)
#' }
#'
#' @importFrom rtracklayer import
#' @importFrom scatterpie geom_scatterpie
#' @importFrom Gviz GeneRegionTrack GenomeAxisTrack DataTrack
#' @importFrom ggplot2 ggplot aes geom_bar theme ylab scale_fill_manual
#'   geom_abline xlab theme_bw expand_limits coord_equal facet_wrap guides
#'   element_text guide_legend
#' @importFrom dplyr filter mutate
#' @importFrom GenomicRanges reduce
#' @importFrom IRanges overlapsAny
#' @importFrom grDevices colorRampPalette dev.off png
#' @importFrom GenomeInfoDb seqnames
#' @importFrom cowplot plot_grid ggdraw draw_image plot_grid get_legend
#' @importFrom S4Vectors %in%
#' @importFrom stats runif
#'
plotGeneSummary <- function(gtf, junctionCovs, useGene, bwFiles,
                            txQuants, txCol = "transcript_id",
                            geneCol = "gene_id", exonCol = "exon_id",
                            cex.title = 0) {
  options(ucscChromosomeNames = FALSE)
  
  genemodels <- rtracklayer::import(gtf)
  idx <- match(c(txCol, geneCol, exonCol), colnames(mcols(genemodels)))
  colnames(mcols(genemodels))[idx] <- c("transcript", "gene", "exon")
  S4Vectors::mcols(genemodels)$symbol <- mcols(genemodels)$transcript
  genemodels <- subset(genemodels, type == "exon")
  
  jl <- junctionCovs %>% dplyr::filter(gene == useGene)
  
  #names(bwFiles) <- ""
  bwConds <- structure(paste0("g", seq_along(bwFiles)), names = names(bwFiles))
  
  ## Gene model track
  if (!("gene_name" %in% colnames(S4Vectors::mcols(genemodels)))) {
    S4Vectors::mcols(genemodels)$gene_name <-
      S4Vectors::mcols(genemodels)[[geneCol]]
  }
  gm <- subset(genemodels, tolower(gene) == tolower(useGene) |
                 tolower(gene_name) == tolower(useGene))
  gm <- subset(gm, gene == gene[1])
  id <- unique(gm$gene_name)
  idshow <- paste0(id, " (", unique(gm$gene), ")")
  show_chr <- unique(seqnames(gm))[1]
  gm <- subset(gm, seqnames == show_chr)
  min_coord <- min(start(gm)) - 0.2*(max(end(gm)) - min(start(gm)))
  max_coord <- max(end(gm)) + 0.05*(max(end(gm)) - min(start(gm)))
  gm$transcript <- factor(gm$transcript, levels = unique(gm$transcript))
  
  ## Other features in the considered region
  gmo <- genemodels[IRanges::overlapsAny(
    genemodels,
    GRanges(seqnames = show_chr,
            ranges = IRanges(start = min_coord,
                             end = max_coord),
            strand = "*"))]
  gmo <- gmo[!(S4Vectors::`%in%`(gmo, gm))]
  gmo <- GenomicRanges::reduce(gmo)
  
  ## Define colors
  muted <- c("#DC050C","#E8601C","#7BAFDE","#1965B0","#B17BA6",
             "#882E72","#F1932D","#F6C141","#F7EE55","#4EB265",
             "#90C987","#CAEDAB","#777777")
  txs <- levels(gm$transcript)
  ncols <- nlevels(gm$transcript)
  cols <- grDevices::colorRampPalette(muted)(ncols)
  names(cols) <- txs
  
  grtr <- Gviz::GeneRegionTrack(gm, showId = TRUE, col = NULL,
                                fill = cols[gm$transcript], name = "",
                                background.title = "transparent",
                                col.title = "black", min.height = 15)
  grtr2 <- Gviz::GeneRegionTrack(gmo, showId = TRUE, col = "black",
                                 fill = "white", name = "", col.title = "black",
                                 showId = FALSE, min.height = 15,
                                 background.title = "transparent")
  
  gtr <- Gviz::GenomeAxisTrack()
  
  tracks <- c(gtr, grtr, grtr2)
  
  multiTracks_rnaseq <- lapply(1:length(bwFiles), function(i) {
    assign(paste0("rnaseqtr", i), 
           Gviz::DataTrack(range = bwFiles[i],
                           name = names(bwFiles)[i], 
                           type = "histogram",
                           chromosome = unique(GenomeInfoDb::seqnames(gm)),
                           col.title = "black",
                           fill = "grey",
                           col = "grey",
                           col.histogram = "grey",
                           fill.histogram = "grey",
                           cex.title = cex.title)
    )})
  tracks <- c(multiTracks_rnaseq, tracks)
  
  rn <- round(1e6*stats::runif(1))
  tmpdir <- tempdir()
  grDevices::png(paste0(tmpdir, "/gviz", rn, ".png"), width = 10.5,
                 height = 5.25, unit = "in", res = 400)
  Gviz::plotTracks(tracks, chromosome = show_chr,
                   from = min_coord, to = max_coord, main = idshow,
                   min.width = 0, min.distance = 0, collapse = FALSE)
  grDevices::dev.off()
  
  tpms <- ggplot2::ggplot(
    txQuants %>% dplyr::filter(gene == useGene) %>%
      dplyr::mutate(
        transcript = factor(transcript,
                            levels = levels(gm$transcript))),
    ggplot2::aes(x = method, y = TPM, fill = transcript)) +
    ggplot2::geom_bar(stat = "identity", position = "fill") +
    ggplot2::ylab("Relative TPM") + ggplot2::xlab("") +
    ggplot2::theme_bw() +
    ggplot2::scale_fill_manual(values = cols, name = "") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5,
                                                       hjust = 1, size = 9),
                   legend.text = ggplot2::element_text(size = 7))
  
  for (tt in txs) {
    jl[[tt]] <- grepl(tt, jl$transcript)
  }
  jl[, txs] <- sweep(jl[, txs], 1, rowSums(jl[, txs]), "/")
  
  jcov <- ggplot2::ggplot() + ggplot2::geom_abline(intercept = 0, slope = 1) +
    scatterpie::geom_scatterpie(ggplot2::aes(x = uniqreads, y = scaledCov,
                                             r = max(scaledCov)/13),
                                cols = txs, data = jl, color = NA) +
    ggplot2::facet_wrap(~ methodscore, nrow = length(unique(jl$methodscore)) %/% 4 + 1) +
    ggplot2::coord_equal(ratio = 1) +
    ggplot2::expand_limits(x = range(c(0, jl$scaledCov, jl$uniqreads)),
                           y = range(c(0, jl$scaledCov, jl$uniqreads))) +
    ggplot2::scale_fill_manual(values = cols, name = "") +
    ggplot2::xlab("Number of uniquely mapped reads spanning junction") +
    ggplot2::ylab("Scaled predicted junction coverage") +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 10),
                   legend.text = ggplot2::element_text(size = 7))
  
  plt <- cowplot::plot_grid(
    cowplot::ggdraw() +
      cowplot::draw_image(paste0(tmpdir, "/gviz", rn, ".png")),
    cowplot::plot_grid(tpms + ggplot2::theme(legend.position = "none"),
                       jcov + ggplot2::theme(legend.position = "none"),
                       nrow = 1, labels = c("B", "C"), rel_widths = c(0.9, 1)),
    cowplot::get_legend(
      tpms + ggplot2::theme(legend.direction = "horizontal",
                            legend.justification = "center",
                            legend.box.just = "bottom") +
        ggplot2::guides(fill = ggplot2::guide_legend(nrow =
                                                       length(txs) %/% 8 + 1))),
    ncol = 1, rel_heights = c(1, 0.7, 0.1),
    labels = c("A", "", ""))
  
  unlink(paste0(tmpdir, "/gviz", rn, ".png"))
  plt
}
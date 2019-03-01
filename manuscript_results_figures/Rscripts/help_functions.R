## Define new panel function for ggpairs, to show both Pearson and Spearman correlation
suppressPackageStartupMessages({
  library(rlang)
})

combinecor <- function(data, mapping, color = I("black"), ...) {
  ## get the x and y data to use the other code
  x <- rlang::eval_tidy(mapping$x, data)
  y <- rlang::eval_tidy(mapping$y, data)
  
  ct1 <- cor(x, y, method = "pearson", use = "pairwise.complete.obs")
  ct2 <- cor(x, y, method = "spearman", use = "pairwise.complete.obs")

  ## print the correlation values
  ggally_text(
    label = paste0("Pearson: \n", signif(ct1, 3), "\nSpearman: \n", signif(ct2, 3)), 
    mapping = mapping,
    xP = 0.5, yP = 0.5, 
    color = color,
    ...
  ) 
}

na_replace <- function(x) {
  x[is.na(x)] <- 0
  x
}


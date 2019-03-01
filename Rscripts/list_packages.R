args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

routdirs <- strsplit(routdirs, ",")[[1]]

print(routdirs)
print(outtxt)

all_packages <- c()

for (rd in routdirs) {
  lf <- list.files(rd, full.names = TRUE)
  for (f in lf) {
    x <- readLines(f)
    idx1 <- which(x == "other attached packages:")
    idx2 <- which(x == "loaded via a namespace (and not attached):")
    if (length(idx1) != 0 & length(idx2) != 0) {
      all_packages <- 
        unique(c(all_packages, 
                 do.call(c, lapply((idx1 + 1):(idx2 - 2), function(i) {
                   grep("\\[", setdiff(setdiff(strsplit(x[i], " ")[[1]], " "), ""), 
                        value = TRUE, invert = TRUE)
                 }))))
    }
  }
}
write.table(sort(all_packages), file = outtxt, 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

## Remove reads that don't have _t or _c in the end of their identifier

print(fastqin)
print(fastqout)

suppressPackageStartupMessages(library(ShortRead))

x <- readFastq(fastqin)
nm <- as.character(x@id)
keep <- grep("_t |_c ", nm)

x <- x[keep]

writeFastq(x, file = fastqout, mode = "w", compress = TRUE)

sessionInfo()
date()

library(DNAcopy)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
input <- args[1]

file_base_name <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(input))
file_base_name <- gsub("_raw.csv", "", file_base_name)

bins <- read_tsv(input, col_types = c("fddd"),
                       col_names = c("chr", "start", "end", "ratio"))

bins <- data.frame(bins)

CNA.object <- CNA(cbind(bins$ratio), bins$chr, bins$start, data.type = "logratio", sampleid = file_base_name)
smoothed.CNA.object <- smooth.CNA(CNA.object)
segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=1, alpha = 0.001)
sdundo.CNA.object <- segment(smoothed.CNA.object,undo.splits="sdundo",undo.SD=3,verbose=1)

segs <- sdundo.CNA.object[["output"]]
#segs <- segment.smoothed.CNA.object[["output"]]
segs <- segs %>% select(-c("ID", "num.mark"))
write_tsv(segs, paste0(file_base_name, "_segments.tsv"))

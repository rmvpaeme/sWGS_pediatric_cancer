library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
#input <- "/Users/rmvpaeme/Repos/RVPCVP2012_sWGS/data/sWGS_healthy/CFD1801932_S12_segments.bed"
file_base_name <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(input))
file_base_name <- gsub("_segments", "", file_base_name)

segs <- read_tsv(input, col_types = c("fdddd"),col_names = TRUE)
segs <- segs %>% filter(!chr %in% c("X", "Y", "23", "24"))
segs <- tibble(segs)


get.cpa <- function(seg){
  cpa <- sum(abs(seg$zscore) * (seg$end - seg$start  + 1)) / nrow(seg) * 1e-8 # per 100 Mb
  return(cpa)
}

get.cpa.modified2 <- function(seg){
  cpa <- sum(abs(seg$ratio) * (seg$end - seg$start + 1)) / 1e8
  return(cpa)
}

line <- data.frame(file = file_base_name, cpa = get.cpa(segs), cpam = get.cpa.modified2(segs))

write_csv(line, paste0(file_base_name, "_cpa.csv"), col_names = FALSE)

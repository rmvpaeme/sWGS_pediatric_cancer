library(tidyverse)

#df <- read_tsv("/Users/rmvpaeme/Dropbox (speleman lab)/Basecamp/CVP/LQB pediatric patients/analysis_RVP/data/IDAT/20201103_850K_MH.txt")
#sample_annotation <- read_csv("/Users/rmvpaeme/Dropbox (speleman lab)/Basecamp/CVP/LQB pediatric patients/analysis_RVP/SNParray_resources/SampleSheet_850K.csv",skip = 4)
#df <- read_tsv("/Users/rmvpaeme/Dropbox (speleman lab)/Basecamp/CVP/LQB pediatric patients/analysis_RVP/data/IDAT/20201103_cytoSNP_MH.txt")
#sample_annotation  <- read_csv("/Users/rmvpaeme/Dropbox (speleman lab)/Basecamp/CVP/LQB pediatric patients/analysis_RVP/SNParray_resources/SampleSheet_cytoSNP.csv",skip = 4)
# df <- read_tsv("/Users/rmvpaeme/Dropbox (speleman lab)/Basecamp/CVP/LQB pediatric patients/analysis_RVP/data/IDAT/20201103_GCT_MH.txt")
# sample_annotation  <- read_csv("/Users/rmvpaeme/Dropbox (speleman lab)/Basecamp/CVP/LQB pediatric patients/analysis_RVP/SNParray_resources/SampleSheet_GCT.csv",skip = 4)
# 

fixed_cols <- 10

df_fixed <- df[,1:10]

nsamples <- (ncol(df) - 10)/4

sample_pre = 11
sample_post <- sample_pre + 4

while (sample_post <= ncol(df)){
  tmp <- df[,sample_pre:sample_post]
  sample_ID <- gsub(".GType", "", colnames(tmp)[1])
  Sample_Name <- sample_annotation %>% filter(Sample_ID == sample_ID) %>% pull(Sample_Name)

  tmp_all <- cbind(df_fixed, tmp)

  write_tsv(tibble(tmp_all), paste0("/Users/rmvpaeme/Dropbox (speleman lab)/Basecamp/CVP/LQB pediatric patients/analysis_RVP/data/IDAT/Extracted/", Sample_Name, "_snpArray.tsv")) 
  
  sample_pre <- sample_post + 1
  sample_post <- sample_post + 5
  }

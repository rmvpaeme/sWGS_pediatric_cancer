library(tidyverse)
library(wesanderson)
library(ggpubr)
library(zoo)
library(readxl)
library(reshape2)
library(ComplexHeatmap)
options(scipen=999)

##########   ##########   ##########   ##########   ##########   ##########   ##########   ##########   ##########   ##########
# bins_top <- read_tsv("/Users/rmvpaeme/test/wiseX.400kb.GRCh38.sWGS/backup/CFD1900265_S2_bins.bed", 
#                      col_types = "fdd?dd")
# 
# segs_top <- read_tsv("/Users/rmvpaeme/test/wiseX.400kb.GRCh38.sWGS/backup/CFD1900265_S2_segments.bed", 
#                      col_types = "fdddd")
# 
# bins_bottom <- read_tsv("/Users/rmvpaeme/test/wiseX.200kb.GRCh38.cfRRBS/snakeCFD1900265-NBL_trimmed_bismark_bt2_bins.bed", 
#                         col_types = "fdd?dd")
# 
# segs_bottom <- read_tsv("/Users/rmvpaeme/test/wiseX.200kb.GRCh38.cfRRBS/snakeCFD1900265-NBL_trimmed_bismark_bt2_segments.bed", 
#                         col_types = "fdddd") ``

setwd("~/Repos/RVPCVP2012_sWGS")
#sample_annotation <- read_excel("sample_annotation_sWGSvsArray.xlsx")
sample_annotation <- read_excel("/Users/rmvpaeme/Dropbox (speleman lab)/Basecamp/CVP/LQB pediatric patients/analysis_RVP/sample_annotation.xlsx")
sample_annotation <- sample_annotation %>% filter(cfDNA_data_available == "TRUE" & tumorDNA_data_available == "TRUE")

chr_order <- c("1", "2", "3", "4", "5", "6", "7", "8",
                               "9", "10", "11", "12", "13", "14", "15",
                               "16", "17", "18", "19", "20", "21", "22", "X", "Y")
               
datadir1 <- "/Users/rmvpaeme/Repos/RVPCVP2012_sWGS/data/sWGS/"
datadir2 <- "/Users/rmvpaeme/Repos/RVPCVP2012_sWGS/data/SNParray/"
sample1 <- "CFD1806833"
sample2 <- "2E71"
binsize <- "200kb"
patientID <- "pt1"
makeCNVcomparison <- function(datadir1, datadir2, sample1, sample2, patientID, binsize, tumor_input){
   
      print(paste0("cfDNA sample: ", sample1))
      print(paste0("tumorDNA sample: ", sample2))
      print(paste0("binsize: ", binsize))
      

      bins_top <- read_tsv(paste0(datadir1, str_subset(dir(datadir1, pattern = sample1), "bins")), col_types = c("fdd?d"),
                                                       col_names = c("chr", "start", "end", "id", "ratio"), skip = 1)
      bins_top <- bins_top %>% filter(!chr %in% c("X", "Y", "23", "24"))
      
      segs_top <- read_tsv(paste0(datadir1, str_subset(dir(datadir1, pattern = sample1), "segments")), col_types = c("fddd"),
                           col_names = c("chr", "start", "end", "ratio"), skip = 1)
      segs_top <- segs_top %>% filter(!chr %in% c("X", "Y", "23", "24")) %>% mutate(Sample = sample1)
      
      if (tumor_input == "array" | tumor_input == "SNParray" ){
         bins_bottom <-  read_tsv(paste0(datadir2,"/",sample2,"_", binsize, ".tsv"), col_types = c("fddd"),
                              col_names = c("chr", "start", "end", "ratio"))
         bins_bottom <- bins_bottom %>% filter(!chr %in% c("X", "Y", "23", "24"))
         bins_bottom$chr <- factor(bins_bottom$chr, levels = chr_order)
         
         segs_bottom <- read_tsv(paste0(datadir2,"/",sample2,"_", binsize, "_segments.tsv"), col_types = c("fddd"),
                              col_names = c("chr", "start", "end", "ratio"), skip = 1)
         
         segs_bottom <- segs_bottom %>% filter(!chr %in% c("X", "Y", "23", "24")) %>% mutate(Sample = sample2)
         segs_bottom$chr <- factor(segs_bottom$chr, levels = chr_order)
      } else if (tumor_input == "sWGS"){
         bins_bottom <-  read_tsv(paste0(datadir2, str_subset(dir(datadir2, pattern = sample2), "bins")), col_types = c("fdd?d"),
                                  col_names = c("chr", "start", "end", "id", "ratio"), skip = 1)
         bins_bottom <- bins_bottom %>% filter(!chr %in% c("X", "Y", "23", "24"))
         bins_bottom$chr <- factor(bins_bottom$chr, levels = chr_order)
         
         segs_bottom <- read_tsv(paste0(datadir2, str_subset(dir(datadir2, pattern = sample2), "segments")), col_types = c("fddd"),
                                 col_names = c("chr", "start", "end", "ratio"), skip = 1)
         segs_bottom <- segs_bottom %>% filter(!chr %in% c("X", "Y", "23", "24")) %>% mutate(Sample = sample2)
         segs_bottom$chr <- factor(segs_bottom$chr, levels = chr_order)
      }
##########   ##########   ##########   ##########   ##########   ##########   ##########   ##########   ##########   ##########   
      # bins_top <- read_tsv("/Users/rmvpaeme/Dropbox (speleman lab)/Basecamp/RVP/RVP1801_cfRRBS/CNV_analysis/bedfiles.wx/CFD1800942_S44_bins.bed", 
      #                      col_types = "fdd?dd")
      # bins_top <- bins_top %>% filter(!chr %in% c("X", "Y", "23", "24"))
      # 
      # segs_top <- read_tsv("/Users/rmvpaeme/Dropbox (speleman lab)/Basecamp/RVP/RVP1801_cfRRBS/CNV_analysis/bedfiles.wx/CFD1800942_S44_aberrations.bed", 
      #                      col_types = "fdddd")
      # segs_top <- segs_top %>% filter(!chr %in% c("X", "Y", "23", "24"))
      # 
# 
#       bins_bottom <- read_tsv("/Users/rmvpaeme/test/wiseX.400kb.GRCh38.sWGS/backup/CFD1806851_S47_bins.bed", col_types = c("fdd?d"),
#                            col_names = c("chr", "start", "end", "id", "ratio"), skip = 1)
#       bins_bottom <- bins_bottom %>% filter(!chr %in% c("X", "Y", "23", "24"))
# 
#       segs_bottom <- read_tsv("/Users/rmvpaeme/test/wiseX.400kb.GRCh38.sWGS/backup/CFD1806851_S47_segments.bed",
#                               col_types = "fddd", col_names = c("chr", "start", "end", "ratio"), skip = 1)
#       segs_bottom <- segs_bottom %>% filter(!chr %in% c("X", "Y", "23", "24"))
# 
#       bins_top <- read_tsv("/Users/rmvpaeme/Repos/RVPCVP2012_sWGS/data/SNParray/2S88_200kb.tsv", col_types = c("fddd"),
#                            col_names = c("chr", "start", "end", "ratio"))
# 
#       bins_top <- bins_top %>% filter(!chr %in% c("X", "Y", "23", "24"))
# 
#            #
#       bins_top$chr <- factor(bins_top$chr, levels = chr_order)
# 
#       segs_top <- read_tsv("/Users/rmvpaeme/Repos/RVPCVP2012_sWGS/data/SNParray/2S88_200kb_segments.tsv", col_types = c("fddd"),
#                            col_names = c("chr", "start", "end", "ratio"), skip = 1)
#       segs_top <- segs_top %>% filter(!chr %in% c("X", "Y", "23", "24"))
#       segs_top$chr <- factor(segs_top$chr, levels = chr_order)
#       
##########   ##########   ##########   ##########   ##########   ##########   ##########   ##########   ##########   ##########  


      color_top <- wes_palette("Cavalcanti1")[1]
      color_bottom <- wes_palette("Cavalcanti1")[4]
      color_abberations <- wes_palette("Cavalcanti1")[5]
      
      max.ratio_bins <- max(c(bins_top$ratio, bins_bottom$ratio), na.rm = TRUE)
      min.ratio_bins <- min(c(bins_top$ratio, bins_bottom$ratio), na.rm = TRUE)
      
      max.ratio_bins.top <- max(c(bins_top$ratio), na.rm = TRUE)
      min.ratio_bins.top <- min(c(bins_top$ratio), na.rm = TRUE)
      
      max.ratio_bins.bottom <- max(c(bins_bottom$ratio), na.rm = TRUE)
      min.ratio_bins.bottom <- min(c(bins_bottom$ratio), na.rm = TRUE)
      
      get.cpa <- function(seg){
         cpa <- sum(abs(seg$zscore) * (seg$end - seg$start  + 1)) / nrow(seg) * 1e-8 # per 1 Mb
         return(cpa)
      }
      
      get.cpa.modified <- function(seg){
         cpa <- sum(abs(seg$ratio) * (seg$end - seg$start  + 1)) / nrow(seg) * 1e-6 # per 1 Mb
         return(cpa)
      }
      
      
      cpa_top <- round(get.cpa.modified(segs_top),4)
      cpa_bottom <- round(get.cpa.modified(segs_bottom),4)
      
      
      ptop <- ggplot(bins_top, aes(x = start, y = ratio)) + 
         theme_bw() + 
         labs(title = paste0(sample1, " CPA: ", cpa_top), y = "log2(ratio)") + 
         geom_point(size = 0.3, col = color_top, alpha = 0.7) +
         lims(y = c(min.ratio_bins.top,max.ratio_bins.top))+
         theme(panel.spacing.x = unit(0, "lines"),
               panel.spacing.y = unit(0, "lines"),
               axis.title.x =  element_blank(),
               axis.text.x =  element_blank(),
               axis.ticks.x =  element_blank(),
               strip.background = element_rect(color = "white", fill = "white"),
               #   panel.border = element_rect(color = "gray", fill = NA, size = 0.3), 
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank()) + 
         geom_hline(yintercept = 0, linetype = "dashed", col = "gray") +
         facet_wrap(~chr, strip.position = "bottom", scales ="free_x", nrow = 1) + 
         geom_segment(data = segs_top, aes(x = start, xend = end, y = ratio , yend = ratio), col = color_abberations) +
         geom_rect(data = segs_top, aes(xmin=start, xmax=end, ymin=0, ymax=ratio), fill = color_abberations, alpha = 0.4) 
      
      
      pbottom <- ggplot(bins_bottom, aes(x = start, y = ratio)) + 
         theme_bw() + 
         labs(title = paste0(sample2, " CPA: ", cpa_bottom), y = "log2(ratio)") +  
         geom_point(size = 0.3, col = color_bottom, alpha = 0.7) +
         lims(y = c(min.ratio_bins,max.ratio_bins))+
         theme(panel.spacing.x = unit(0, "lines"),
               panel.spacing.y = unit(0, "lines"),
               axis.title.x =  element_blank(),
               axis.text.x =  element_blank(),
               axis.ticks.x =  element_blank(),
               strip.background = element_rect(color = "white", fill = "white"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank()) + 
         geom_hline(yintercept = 0, linetype = "dashed", col = "gray") +
         facet_wrap(~chr, strip.position = "bottom", scales ="free_x", nrow = 1) +
         geom_segment(data = segs_bottom, aes(x = start, xend = end, y = ratio , yend = ratio), col = color_abberations) +
         geom_rect(data = segs_bottom, aes(xmin=start, xmax=end, ymin=0, ymax=ratio), fill = color_abberations, alpha = 0.4) 
      
      
      corplot_df_bins <- full_join(bins_top, bins_bottom, by = c("chr", "start", "end"))
      #corplot_df <- corplot_df %>% filter(!is.na(ratio.x) & !is.na(ratio.y))
      
      corplot_df_bins <- corplot_df_bins %>%
         mutate(rM.top=rollapply(ratio.x,50, FUN=function(x) mean(x, na.rm=TRUE), fill=NA, align="right")) %>%
         mutate(rM.bottom=rollapply(ratio.y,50, FUN=function(x) mean(x, na.rm=TRUE), fill=NA, align="right")) 
      
      rollmean <- ggplot(tibble(corplot_df_bins), aes(x = start)) + 
         theme_bw() + 
         labs(y = "log2(ratio)", x = "chromosomes") +
         geom_line(aes(y=rM.bottom), col = color_bottom, alpha = 1, size = 1) +
         geom_line(aes(y=rM.top), col = color_top, alpha = 0.7, size = 1) +
         theme(panel.spacing.x = unit(0, "lines"),
               panel.spacing.y = unit(0, "lines"),
               #  axis.title.x =  element_blank(),
               axis.text.x =  element_blank(),
               axis.ticks.x =  element_blank(),
               strip.background = element_rect(color = "white", fill = "white"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank()) + 
         geom_hline(yintercept = 0, linetype = "dashed", col = "gray") +
         facet_wrap(~chr, strip.position = "bottom", scales ="free_x", nrow = 1) + lims(y = c(-NA,NA))
      
      max.ratio <- max(c(corplot_df_bins$ratio.x, corplot_df_bins$ratio.y), na.rm = TRUE)
      min.ratio <- min(c(corplot_df_bins$ratio.x, corplot_df_bins$ratio.y), na.rm = TRUE)
      corplot <- ggplot(corplot_df_bins, aes(x = ratio.x, y = ratio.y)) + 
         geom_point(size = 0.3, col = wes_palette("Royal1")[1]) + 
         theme_bw() +
         geom_abline(slope = 1, intercept = 0) +
         stat_cor(method = "pearson", aes(label = ..r.label..)) +
         geom_smooth(method = "lm", col = color_abberations, linetype = "dashed") +
         lims(x = c(min.ratio , max.ratio), y = c(min.ratio, max.ratio))
      
      R <- cor(corplot_df_bins$ratio.x, corplot_df_bins$ratio.y, use = "complete.obs", method = "pearson")
      
      CNVplot <- ggarrange(ptop, pbottom, nrow = 2, ncol = 1, align = "v", heights = c(1, 1))
      comparisons <- ggarrange(rollmean, corplot, nrow = 1, ncol = 2, widths = c(0.6, 0.4))
      
      arrangeplot <- ggarrange(CNVplot, comparisons, nrow = 2, heights = c(0.75, 0.25), widths = c(1, 0.90), align = "hv")
      ggsave(paste0("/Users/rmvpaeme/Repos/RVPCVP2012_sWGS/plots/", sample1, "_", sample2, "_", binsize, ".png"), plot = arrangeplot, width = 12, height = 10, dpi = 300)
      
      df_comparison <- data.frame(PatientID = c(patientID), pearsonR = c(R), cpa_cfDNA = c(cpa_top), cpa_tumorDNA = c(cpa_bottom))
      return(df_comparison)
}


cfDNAvsTissue <- data.frame()
sample_annotation <- sample_annotation  %>% filter(SampleOrigin == "MH")
for (row in 1:nrow(sample_annotation)){
   sample1 <- sample_annotation[row,]$CFD_ID
   sample2 <- sample_annotation[row,]$tumorDNA_ID
   patientID <-  sample_annotation[row,]$PatientID
   datadir1 <- "/Users/rmvpaeme/Repos/RVPCVP2012_sWGS/data/sWGS/"
   if (sample_annotation[row,]$tumorDNA_assay == "array"){
      input_tumor = "array"
      datadir2 <-  "/Users/rmvpaeme/Repos/RVPCVP2012_sWGS/data/array/"
   } else if (sample_annotation[row,]$tumorDNA_assay == "SNParray"){
      input_tumor = "SNParray"
      datadir2 <-  "/Users/rmvpaeme/Repos/RVPCVP2012_sWGS/data/SNParray/"
   } else if (sample_annotation[row,]$tumorDNA_assay == "sWGS"){
      input_tumor = "sWGS"
      datadir2 <-  "/Users/rmvpaeme/Repos/RVPCVP2012_sWGS/data/sWGS/"
   }
   if (length(dir(datadir2, pattern = sample2) > 0  | length(dir(datadir1, pattern = sample1) > 0))) {
      tmp <- makeCNVcomparison(datadir1, datadir2, sample1, sample2, patientID, "200kb", input_tumor)
      cfDNAvsTissue <- rbind(cfDNAvsTissue, tmp)
   }
}

cfDNAvsTissue <- merge(sample_annotation, cfDNAvsTissue, by = "PatientID")
ggplot(cfDNAvsTissue, aes(y = as.numeric(cpa_cfDNA), x = cpa_tumorDNA)) + geom_point()


read_segments <- function(datadir){
   files<-list.files(c(datadir),recursive=TRUE)
   files<-files[grep("segments_per", files)]
   tmp_heatmap <- data.frame()
   for(i in 1:length(files)){
         name_sample <- gsub("_segments_per_[0-9]+kb.tsv", "", files[i])
         print(name_sample)
         tmp <-  read_tsv(paste0(datadir,files[i]), col_types = c("fddd"),
                          col_names = c("chr", "start", "end", "ratio"), skip = 1)
      
         colnames(tmp)<- c("chr", "start", "end", "ratio")
         tmp$SampleID <- name_sample
         if (nrow(tmp_heatmap) == 0){
            tmp_heatmap <- tmp
         } else {
         tmp_heatmap <- rbind(tmp_heatmap, tmp)
         }
   }
   return(tmp_heatmap)
}
df_heatmap_init <- data.frame()
df_heatmap_init <- rbind(df_heatmap_init, read_segments(datadir1))
df_heatmap_init <- rbind(df_heatmap_init, read_segments(datadir2))

df_heatmap <- df_heatmap_init
df_heatmap$bin <- paste0(as.character(df_heatmap$chr), ":", as.character(df_heatmap$start), "-", as.character(df_heatmap$end))

df_heatmap <- df_heatmap %>% select(-c(chr, start, end)) %>% spread(key = c(bin), value = ratio)
df_heatmap$SampleID <- gsub("_.*", "", df_heatmap$SampleID)
#df_heatmap <- df_heatmap %>% distinct(SampleID, .keep_all = TRUE)

colnames_hm <- paste0(gsub(":.*", "", colnames(df_heatmap[,2:ncol(df_heatmap)])))
colnames_hm <- factor(colnames_hm, levels = chr_order)
rownames_hm <- df_heatmap$SampleID

tumorSamples <- sample_annotation %>% select(UniqueID, PatientID, tumorDNA_ID, TumorType,tumorDNA_assay) %>% filter(tumorDNA_ID %in% rownames_hm)
colnames(tumorSamples) <- c("UniqueID", "PatientID",  "SampleID", "TumorType", "assay")
tumorSamples$biomaterial <- "tumor DNA"

cfDNASamples <-  sample_annotation %>% select(UniqueID, PatientID, CFD_ID, TumorType,tumorDNA_assay) %>% filter(CFD_ID %in% rownames_hm)
colnames(cfDNASamples) <- c("UniqueID", "PatientID",  "SampleID", "TumorType", "assay")
cfDNASamples$assay <- "sWGS"
cfDNASamples$biomaterial <- "cfDNA"

sampleTypes <- bind_rows(tumorSamples, cfDNASamples)
sampleTypes <- sampleTypes %>% dplyr::arrange(PatientID)

df_heatmap <- inner_join(sampleTypes, df_heatmap)

df_heatmap <- df_heatmap %>% dplyr::arrange(TumorType, patientID)
#    cbind(sample_annotation$sWGS, sample_annotation$TumorType, rep("sWGS", nrow(sample_annotation))),
#    cbind(sample_annotation$array, sample_annotation$TumorType, rep("array", nrow(sample_annotation)))
# ))
# colnames(sampleTypes) <- c("SampleID", "tumor", "type")


ha = HeatmapAnnotation(
   platform = df_heatmap$assay,
   biomaterial = df_heatmap$biomaterial,
   tumor = df_heatmap$TumorType,
   annotation_name_side = "bottom", which = "row", show_legend = TRUE,
   width = unit(0.7, "cm"),
   show_annotation_name = TRUE,
   annotation_legend_param = 
      list(
         platform = list(
            title = "platform",
            at = c("SNParray", "sWGS"),
            labels = c("SNParray", "sWGS"),
            col = c("pink", "navyblue")
         )
      ))


ht <- Heatmap(as.matrix(df_heatmap[,7:ncol(df_heatmap)]), 
        name = "log2ratio", column_split = colnames_hm, row_split = df_heatmap$PatientID,
        row_title = NULL,
        row_gap = unit(0.2, "mm"),
        cluster_columns = FALSE, 
        cluster_rows = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE, 
        row_labels = rownames_hm, 
        right_annotation = ha,
        column_title_gp = gpar(fontsize = 9),
        column_gap = unit(1, "mm"))
draw(ht,heatmap_legend_side = "bottom")

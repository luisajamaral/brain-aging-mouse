
#Great enrichment 

library(rGREAT)
library(GenomicRanges)

setwd("/mnt/tscc/lamaral/projects/brain_aging_mouse/analysis/snapATAC/all66/")


getGREATTableBG <- function(file,cl) {
  
  
  bedfile = read.csv(file, sep = "\t",header = F)
  if (length(grep("rand", bedfile$V1))>0) {bedfile = bedfile[-grep("rand", bedfile$V1),]}
  if (length(grep("Un", bedfile$V1))>0) {bedfile = bedfile[-grep("Un", bedfile$V1),]}
  if (length(grep("M", bedfile$V1))>0) {bedfile = bedfile[-grep("M", bedfile$V1),]}
  
  if(nrow(bedfile))
  
  bg = read.csv(paste(cl,"_cpm.txt",sep = ""),header = T)
  bg_coors = bg$V1
  chr = sapply(strsplit(as.character(bg$X), ":"), "[[", 1)
  coors = sapply(strsplit(as.character(bg$X), ":"), "[[", 2)
  start = sapply(strsplit(as.character(coors), "-"), "[[", 1)
  end = sapply(strsplit(as.character(coors), "-"), "[[", 2)
  bg = data.frame(chr,start ,end)
  if (length(grep("rand", bg$chr))>0) {bg = bg[-grep("rand", bg$chr),]}
  if (length(grep("Un", bg$chr))>0) {bg = bg[-grep("Un", bg$chr),]}
  if (length(grep("M", bg$chr))>0) {bg = bg[-grep("M", bg$chr),]}
  bg$start=as.numeric(as.character(bg$start))
  bg$end=as.numeric(as.character(bg$end))
  
  
  
  
  which(!(paste(bedfile$V1,bedfile$V2,bedfile$V3 ) %in% paste(bg$chr,bg$start,bg$end)))
  
    #show(head(bg))
  job = submitGreatJob(bedfile,species = "mm10", request_interval = 1, bg = bg)
  e = getEnrichmentTables(job)
  e = e$`GO Biological Process`
  e =e[order(e$Hyper_Adjp_BH),]
  show(head(e))
  return(e)
}


dat_up_bg_filt = list()
dat_up_bg = list()
files = list.files(".",pattern = "_cluster_peaks.txt")
for (file in files) {
  cl = strsplit(file, "_")[[1]][1]
  color = strsplit(file, "_")[[1]][2]
  cat(file,"\n")
  dat_up_bg[[cl]] = getGREATTableBG(file,cl)
  e = dat_up_bg[[cl]]
  #dat_up_bg_filt[[i]] = e[which(e$Binom_Observed_Region_Hits>5 & e$Binom_Fold_Enrichment>=1 & e$Hyper_Adjp_BH < 0.05),]
  dat_up_bg_filt[[cl]] = e[which(e$Hyper_Total_Regions>5 & e$Hyper_Fold_Enrichment>=1.5 & e$Hyper_Foreground_Gene_Hits>= 3 & e$Hyper_Adjp_BH < 0.05),]
  
  if(nrow(dat_up_bg_filt[[cl]])>0) {
    write.table(dat_up_bg_filt[[cl]], file = paste("great/",cl,"_",color,"_GREAT_05.txt", sep = ""), quote = F, sep = "\t")
  }
}


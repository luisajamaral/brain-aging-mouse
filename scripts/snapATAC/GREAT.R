#Great enrichment 

library(rGREAT)
library(GenomicRanges)

#setwd("/mnt/tscc/lamaral/projects/brain_aging_mouse/analysis/snapATAC/9mo+18mo/age_diff_edgeR.snap_celltype/")
setwd("/mnt/tscc/lamaral/projects/brain_aging_mouse/analysis/snapATAC/all66/age_diff_edgeR.snap/")

num_cell_types = 21

getGREATTable <- function(file) {
  bedfile = read.csv(file, sep = "\t",header = F)
  job = submitGreatJob(bedfile,species = "mm10", request_interval = 1)
  e = getEnrichmentTables(job)
  e = e$`GO Biological Process`
  e =e[order(e$Hyper_Adjp_BH),]
  e= e[which(e$Hyper_Total_Genes>5 & e$Hyper_Fold_Enrichment>=1.5 & e$Hyper_Adjp_BH < 0.1),]
  head(e)
  return(e)
}

getGREATTableBG <- function(file,bg) {
  bedfile = read.csv(file, sep = "\t",header = F)
  
  if (length(grep("rand", bedfile$V1))>0) {bedfile = bedfile[-grep("rand", bedfile$V1),]}
  if (length(grep("Un", bedfile$V1))>0) {bedfile = bedfile[-grep("Un", bedfile$V1),]}
  if (length(grep("M", bedfile$V1))>0) {bedfile = bedfile[-grep("M", bedfile$V1),]}
  
  
  #show(head(bg))
  job = submitGreatJob(bedfile,species = "mm10", request_interval = 1, bg = bg)
  e = getEnrichmentTables(job)
  e = e$`GO Biological Process`
  e =e[order(e$Hyper_Adjp_BH),]
  show(head(e))
  return(e)
}

bg_file = "~/projects/brain_aging_mouse/data/snATAC/peaks/9mo+18mo_summits_ext1k.bed"
bg = read.csv(bg_file, sep = "\t",header = F)
bg = bg[-grep("rand", bg$V1),]
bg = bg[-grep("Un", bg$V1),]
bg = bg[-grep("M", bg$V1),]

dat_up_bg_filt = list()
dat_up_bg = list()
for (i in 1:num_cell_types) {
  cat(i,"\n")
  file = paste("up_bed/",i,".edger.up.bed", sep = "")
  dat_up_bg[[i]] = getGREATTableBG(file,bg)
  e = dat_up_bg[[i]]
  dat_up_bg_filt[[i]] = e[which(e$Hyper_Total_Regions>5 & e$Hyper_Fold_Enrichment>=1.5 & e$Hyper_Foreground_Gene_Hits>= 3 & e$Hyper_Adjp_BH < 0.1),]
  if(nrow(dat_up_bg_filt[[i]])>0) {
    write.table(dat_up_bg_filt[[i]], file = paste("great/",i,"_up_GREAT_01.txt", sep = ""), quote = F, sep = "\t")
  }
}

dat_down_bg = list()
dat_down_bg_filt = list()
for (i in 1:num_cell_types) {
  cat(i,"\n")
  file = paste("down_bed/",i,".edger.down.bed", sep = "")
  dat_down_bg[[i]] = getGREATTableBG(file,bg)
  e = dat_down_bg[[i]]
  dat_down_bg_filt[[i]] = e[which(e$Hyper_Total_Regions>5 & e$Hyper_Fold_Enrichment>=1.5 & e$Hyper_Foreground_Gene_Hits>= 3 & e$Hyper_Adjp_BH < 0.1),]
  if(nrow(dat_down_bg_filt[[i]])>0) {
    write.table(dat_down_bg_filt[[i]], file = paste("great/",i,"_down_GREAT_01.txt", sep = ""), quote = F, sep = "\t")
  }
}





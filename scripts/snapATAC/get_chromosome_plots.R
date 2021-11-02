library(chromPlot)
data("mm10_gap")
mm10_gap[384,] = c("chrY", 0, 100000, "telomere")

setwd("~/projects/brain_aging_mouse/analysis/snapATAC/all66/age_diff_edgeR.snap_final_clusters")
files = list.files(".", "*9mo_vs_18mo.edger.txt")
pdf("chromosome_location_DARs_plot_final_9mo_vs_18mo_fdr.pdf",
    height = 5.5,
    width = 6.5)
for (f in files) {
  diff = fread(f, sep = "\t")
  if (nrow(diff) == 0) {
    next()
  }
  f = gsub("_9mo_vs_18mo.edger",  "", f)
  f = gsub(".txt", "", f)

  head(diff)
  coors = as.character(diff$V1)
  chr = sapply(strsplit(coors, ":"), "[[", 1)
  coors = sapply(strsplit(coors, ":"), "[[", 2)
  start = sapply(strsplit(coors, "-"), "[[", 1)
  end = sapply(strsplit(coors, "-"), "[[", 2)
  
  diff$Chrom = chr
  diff$Start = start
  diff$End = end
  diff$sig = -log10(diff$fdr)
  diff$sig[which(diff$sig>300)] = 300
  #if (length(which(diff$fdr < 0.05))>100) {
  diff = diff[which(diff$fdr < 0.05), ]
#  } else {
#    diff = diff[1:100,]
#  }
  if(nrow(diff)==0){
    next()
  }
  diff$Name = ""
  up = diff[which(diff$logFC < 0), ]
  up = data.frame(up$Chrom, up$Start, up$End, up$sig, up$Name)
  colnames(up) = c("Chrom", "Start", "End", "sig", "Name")
  down = diff[which(diff$logFC > 0), ]
  down = data.frame(down$Chrom, down$Start, down$End, down$sig, down$Name)
  colnames(down) = c("Chrom", "Start", "End", "sig", "Name")

  up$Chrom = as.character(up$Chrom)
  up$Start = as.integer(as.character(up$Start))
  up$End = as.integer(as.character(up$End))
  
  down$Chrom = as.character(down$Chrom)
  down$Start = as.integer(as.character(down$Start))
  down$End = as.integer(as.character(down$End))

  if(nrow(up)<4) {
    up = rbind(up,up)
    up = rbind(up,up)
  }
  
  if(nrow(down)<4) {
    down = rbind(down,down)
  }
  
  if (nrow(up) > 1 & nrow(down) > 1) {
    chromPlot(
      gaps = mm10_gap,
      bands = mm10_cytoBandIdeo,
      stat = up,
      stat2 = down,
      statCol = "sig",
      statCol2 = "sig",
      statName = "up",
      statName2 = "down",
      colStat = "red",
      colStat2 = "blue",
      statTyp  = "p",
      scex = 2,
      spty = 20,
      statThreshold = 1.30103,
      statThreshold2 = 1.30103,
     # chr=c("X"),
      bin = 1e10,
      cex = 1.25,yAxis = F,
      statSumm = "none",
      stack = FALSE,
      title = f
    )
  } else if (nrow(up) > 1) {
    chromPlot(
      gaps = mm10_gap,
      bands = mm10_cytoBandIdeo,
      stat = up,
      statCol = "sig",
      statName = "up",
      colStat = "red",
      statTyp  = "p",
      scex = 2,
      spty = 20,
      statThreshold = 1.30103,
      statThreshold2 = 1.30103,
      # chr=c("chr13"),
      bin = 1e6,
      cex = 0.7,
      statSumm = "none",
      stack = FALSE,
      title = f
    )
  } else if (nrow(down) > 1) {
    chromPlot(
      gaps = mm10_gap,
      bands = mm10_cytoBandIdeo,
      stat = down,
      statCol = "sig",
      statCol2 = "sig",
      statName = "down",
      colStat = "blue",
      statTyp  = "p",
      scex = 2,
      spty = 20,
      statThreshold = 1.30103,
      statThreshold2 = 1.30103,
      chrSide = c(-1, -1, -1, -1, -1, -1, 1, -1),
     # chr=c("Y"),yAxis = T,
      bin = 1e6,
      cex = 0.7,
      statSumm = "none",
      stack = FALSE,
      title = f
    )
}
}
dev.off()

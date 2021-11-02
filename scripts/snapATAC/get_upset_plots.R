
library(UpSetR)

files = list.files("age_diff_edgeR.snap",pattern = ".txt",full.names = T)
cpmfiles = list.files(".",pattern="cpm.txt")

for (f in files){
  tmp = read.csv(f)
  sig = tmp[which(tmp$fdr < 0.05), ]
  #  sig = tmp[which(tmp$PValue< quantile(tmp$PValue,0.01)),]
  chr = sub("(chr.*):(.*)-(.*)", "\\1", sig$X)
  start = sub("(chr.*):(.*)-(.*)", "\\2", sig$X)
  end = sub("(chr.*):(.*)-(.*)", "\\3", sig$X)
  out = data.frame(chr, start, end, sig$logFC, as.integer(-log10(sig$PValue)))
  #  write.table(out, paste0("../age_diff_edgeR.snap_celltype/",cl,".both.bed"),row.names=F,col.names=F,quote=F,sep="\t")
  outF =gsub(".txt", "", f)
  write.table(
    subset(out, sig.logFC > 0),
    paste0(outF,".up.bed"),
    row.names = F,
    col.names = F,
    quote = F,
    sep = "\t"
  )
  write.table(
    subset(out, sig.logFC < 0),
    paste0(outF, ".down.bed"),
    row.names = F,
    col.names = F,
    quote = F,
    sep = "\t"
  )
  
}


pdf("Down_in_aging_upset_plots.pdf",height = 5, width = 5)
for(cl in unique(clusters)) {
  files = list.files("age_diff_edgeR.snap/up_bed",pattern = paste("^",cl,"_",sep = ""),full.names = T)
  peaks = list()
  for(f in files) {
    n = gsub("age_diff_edgeR.snap/up_bed/","", f)
    n = gsub(".edger.up.bed","", n)
    if(length(grep("8wk_vs_9mo",f))>0) {
      f = gsub("up","down",f)
    }
    tmp = read.table(f)
    tpeaks = paste(tmp$V1,":",tmp$V2,"-",tmp$V3, sep = "")
    peaks[[n]] = tpeaks
  }
  print(upset(fromList(peaks),text.scale = 1.3,mainbar.y.label = "# DAR (fdr<0.05)"))
}
dev.off()

clusters = unique(cluster)[-c(11,13,19,17,18,15)]

pdf("Up_in_aging_upset_plots.pdf", height  = 5, width =5)
for(cl in unique(clusters)) {
  files = list.files("age_diff_edgeR.snap/down_bed",pattern = paste("^",cl,"_",sep = ""),full.names = T)
  peaks = list()
  for(f in files) {
    n = gsub("age_diff_edgeR.snap/down_bed/","", f)
    n = gsub(".edger.down.bed","", n)
    if(length(grep("8wk_vs_9mo",f))>0) {
      f = gsub("down","up",f)
    }
    tmp = read.table(f)
    tpeaks = paste(tmp$V1,":",tmp$V2,"-",tmp$V3, sep = "")
    peaks[[n]] = tpeaks
  }
  print(upset(fromList(peaks),text.scale = 1.3,mainbar.y.label = "# DAR (fdr<0.05)"))
}
dev.off()


clusters = unique(cluster)[-c(11,13,19,17,18,15,14,16)]

for(cl in unique(clusters)[16]) {
  cat(cl,"\n")
  f1 = list.files("age_diff_edgeR.snap/down_bed",pattern = paste("^",cl,"_",sep = ""),full.names = T)
  f2 = list.files("age_diff_edgeR.snap/up_bed",pattern = paste("^",cl,"_",sep = ""),full.names = T)
  files = c(f1,f2)
  peaks = c()
  for(f in files) {
    tmp = read.table(f)
    tpeaks = paste(tmp$V1,":",tmp$V2,"-",tmp$V3, sep = "")
    peaks = c(peaks,tpeaks)
  }
  peaks = unique(peaks)
  sig = 1
  while (length(peaks)>15000){
    sig = sig+1
    peaks = c()
    for(f in files) {
      tmp = read.table(f)
      if (length(grep("9mo_vs_18mo", f))>0) {
        tmp = tmp[which(tmp$V5>1),]
      } else {
        tmp = tmp[which(tmp$V5>sig),]
      }
      tpeaks = paste(tmp$V1,":",tmp$V2,"-",tmp$V3, sep = "")
      peaks = c(peaks,tpeaks)
      peaks=unique(peaks)
    }
  }
  cat(sig, "\n")
  cpm = read.csv(paste(cl,"_cpm.txt", sep = ""), row.names = 1,check.names = F)
  peaks = peaks[which(peaks%in%rownames(cpm))]
  if (length(peaks) > 2) {
    o = cpm[peaks,]
    hr <- hclust(as.dist(1-cor(t(o))),method = "ward.D2")
    mycl <- cutree(hr, h = 18)
    c25 <- c(
      "dodgerblue2", "#E31A1C", # red
      "green4",
      "#6A3D9A", # purple
      "#FF7F00", # orange
      "gold1",
      "maroon", "orchid1", "deeppink1","green1", "steelblue4",
      "darkturquoise",  "yellow4", "yellow3",
      "darkorange4", "brown", "blue1"
    )
    
    mycolhc <- c25[1:length(unique(mycl))]
    mycolhc <- mycolhc[as.vector(mycl)]
    o = o[,order(colnames(o))]
    oder = c(colnames(o)[grep("8wk", colnames(o))],colnames(o)[grep("9mo", colnames(o))],colnames(o)[grep("18mo", colnames(o))])
    o = o[,oder]
    col = substr(colnames(o), 1,3)
    #viridis()
    col[which(col=="8wk")] = heat.colors(n = 3)[3]
    col[which(col=="9mo")] = heat.colors(n = 3)[2]
    col[which(col=="18m")] = heat.colors(n = 3)[1]
    
    main = paste("cluster ",cl,"\n", "fdr<10^-", sig, " (",length(peaks)," peaks)",sep = "")
    my_palette =colorpanel(100, "darkblue", "white", "red")
    pdf(paste(cl, "groupedheat.pdf"), height = 7,width = 10)
    heatmap.2(
      data.matrix(o),
      cexCol = 1,
      margins = c(15, 10),
      Rowv = as.dendrogram(hr),
      Colv = F,dendrogram = "row",
      col = my_palette,
      scale = "row",
      trace = "none",
      RowSideColors = mycolhc,
      ColSideColors = col,labRow = F,
      main = main
    )
    legend("topright", title = "age",legend=c("8wk","9mo", "18mo"), 
           fill=c(heat.colors(n = 3)[3],heat.colors(n = 3)[2],heat.colors(n = 3)[1]), cex=0.8, box.lty=0)
    dev.off()
    peaks = rownames(o)
    chr = sapply(strsplit(as.character(peaks), ":"), "[[", 1)
    coors = sapply(strsplit(as.character(peaks), ":"), "[[", 2)
    start = sapply(strsplit(as.character(coors), "-"), "[[", 1)
    end = sapply(strsplit(as.character(coors), "-"), "[[", 2)
    pbed = data.frame(chr,start ,end)
    
    write.table(pbed[which(mycl == 1),],file = paste(cl,"_Blue_cluster_peaks.txt", sep = ""), sep = "\t", col.names = F,row.names = F, quote = F)
    if(2%in%mycl) {
      write.table(pbed[which(mycl == 2),],file = paste(cl,"_Red_cluster_peaks.txt", sep = ""), sep = "\t", col.names = F,row.names = F, quote = F)
    }
    if(3%in%mycl) {
      write.table(pbed[which(mycl == 3),],file = paste(cl,"_Green_cluster_peaks.txt", sep = ""), sep = "\t", col.names = F,row.names = F, quote = F)
    }
    if (4%in%mycl) {
      write.table(pbed[which(mycl == 4),],file = paste(cl,"_Purple_cluster_peaks.txt", sep = ""), sep = "\t", col.names = F,row.names = F, quote = F)
    }
  }
}

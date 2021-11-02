setwd("../../analysis/snapATAC/9mo+18mo/")
x.sp = load("all.cluster.RData")
region = gsub("_2$","",x.sp@metaData$sample)
region = gsub("_1$","",region)
region = gsub("18mo_","",region)
region = gsub("9mo_","",region)
x.sp@metaData$region = region

tab = data.frame(sample=x.sp@sample,cluster=x.sp@cluster,region = x.sp@metaData$region)
summ = plyr::count(tab[,1])
tab2 = plyr::count(tab)
tab2$total = summ$freq[match(tab2$sample,summ$x)]
tab2$frac = tab2$freq/tab2$total

pdf("celltype_frac_per_region.pdf",height=5,width=10)
for(i in unique(x.sp@metaData$region)) {
  p = ggplot(tab2[which(tab2$region==i),]) + 
    geom_col(aes(x=cluster,y=frac,fill=sample),position="dodge") +
    theme_bw() +ggtitle(i)
  print(p)
}
dev.off()




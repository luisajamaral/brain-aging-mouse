library(ggplot2)
tissue = "all66"
RData = "all.cluster.meta.final.RData"
setwd(paste0("../../analysis/snapATAC/", tissue))

load(RData)

sample = paste(x.sp@metaData$stage,"_",x.sp@metaData$region,"_",x.sp@metaData$replicate, sep = "")
sample = gsub(" ",".", sample)
region = x.sp@metaData$region

tab = data.frame(sample=sample,cluster=x.sp@metaData$final_clusters,region = region,rep = x.sp@metaData$replicate,age = x.sp@metaData$stage,
                 age_rep = paste(x.sp@metaData$stage, "_",x.sp@metaData$replicate,sep = ""))
summ = plyr::count(tab[,1])
tab2 = plyr::count(tab)
tab2$total = summ$freq[match(tab2$sample,summ$x)]
tab2$frac = tab2$freq/tab2$total
tab2$age = factor(tab2$age, levels = c("8wk","9mo","18mo"),ordered = T)
pdf("celltype_frac_per_region_fi.pdf",height=5,width=10)
for(i in unique(region)) {
    p = ggplot(tab2[which(tab2$region==i),]) + 
    geom_col(aes(x=cluster,y=frac,fill = age,color=rep),position="dodge") +
    theme_bw() +ggtitle(i) + theme(axis.text.x = element_text(angle = 90))
  print(p)
}
dev.off()


pdf("celltype_frac_per_cluster_fi.pdf",height=5,width=10)
for(i in unique(tab2$cluster)) {
  p = ggplot(tab2[which(tab2$cluster == i),]) + 
    geom_col(aes(x=region,y=frac,fill = age,color=rep),position="dodge") +
    theme_bw() +ggtitle(i) + theme(axis.text.x = element_text(angle = 90))
  print(p)
}
dev.off()
m = x.sp@metaData
m$stage = factor(m$stage, levels = c("8wk", "9mo","18mo"),ordered = T)
m = m[-which(is.na(m$final_clusters)),]
m$final_clusters = factor(m$final_clusters, levels = unique(m$final_clusters))
pdf("num_cells_per_cluster_fi.pdf",height = 5, width =13)
ggplot(m,aes(x=final_clusters,fill = stage)) + 
  geom_bar(stat = "count",position = "dodge") +theme(axis.text.x = element_text(angle = 90))+ geom_hline(yintercept=2000)
  theme_bw() #+ggtitle(i) + theme(axis.text.x = element_text(angle = 90))
dev.off()

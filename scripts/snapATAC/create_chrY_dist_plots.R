setwd("~/projects/brain_aging_mouse/analysis/snapATAC/all66/")
RData = "all.cluster.meta.final.pmat.RData"
load(RData)
peak = x.sp@peak$name
x.sp@metaData$final_clusters[which(x.sp@metaData$final_clusters %in% c("D1MSN","D2MSN"))] = "MSN"
cluster = x.sp@metaData$final_clusters
cluster = gsub("/","-",cluster)
sample = x.sp@metaData$sample_name
stage = x.sp@metaData$stage
tissue = x.sp@metaData$region
tissue = gsub(" ", ".", tissue)
cluster_sample = paste(cluster, sample, sep = ".")
age_cluster_tissue_rep = paste(cluster,stage,tissue, x.sp@metaData$replicate, sep = "_")
age_tissue_rep = paste(stage,tissue, x.sp@metaData$replicate, sep = "_")
cluster_tissue = paste(cluster,tissue,sep = "_")
pmat = x.sp@pmat 
peaks = x.sp@peak$name
pmat = pmat/rowSums(pmat)
chr = sapply(strsplit(as.character(peaks), ":"), "[[", 1)
pmat = pmat[, which(chr == 'chrY')]
rm(x.sp)
unique(cluster)
pdf("chrY_count_dist.pdf",height =4, width =4)
for (cl in unique(cluster)[-29]) { 
mat = pmat[which(cluster == cl), ]
c = rowSums(mat)
c_stage = stage[which(cluster == cl)]

c = data.frame(c_stage, c)
colnames(c) = c("stage","c")
c$stage = factor(c$stage, levels = c("18mo","9mo","8wk"))
# g = ggplot(c, aes(x=c, color=stage)) +
#   geom_histogram(fill="white", alpha=0.1, position="identity",bins =50) + ggtitle(paste("chrY\n", cl))
# print(g)

cdat <- ddply(c, "stage", summarise, rating.mean=mean(c))
g = ggplot(c, aes(x=c, colour=stage)) +
  geom_density() +
  geom_vline(data=cdat, aes(xintercept=rating.mean,  colour=stage),
             linetype="dashed", size=.5)  + xlab(label = "chrY normalized count sum")+ 
  ggtitle(paste(cl),subtitle = "(chrY)")
print(g)

p<-ggplot(c, aes(x=stage, y=c, color=stage)) +
  geom_violin(trim=FALSE) + ggtitle(paste(cl),subtitle = "(chrY)")+
  geom_boxplot(width=0.1,notch = T)+ylab(label = "chrY normalized count sum")
print(p)

}
dev.off()

c = rowSums(pmat)
c = cbind(c,paste(stage), cluster)
colnames(c) = c("c","stage","cluster")
c = as.data.frame(c)

head(c)
stage = factor(stage, levels = c("18mo","9mo","8wk"))
c$c = as.numeric(c$c)
c = c[-which(is.na(c$cluster)),]
pdf("all_clusters_chrY_boxplot.pdf",height = 4, width = 10)
ggplot(c, aes(x=cluster, y=c, fill=stage)) +
  geom_boxplot(notch = T) +  ylim(c(0,0.002))+ theme(axis.text.x = element_text(angle = 90))
dev.off()

ggplot(c, aes(x=cluster, y=c, fill=stage)) +
  geom_violin(trim = T,) + ylim(c(0,0.002))+theme(axis.text.x = element_text(angle = 90))


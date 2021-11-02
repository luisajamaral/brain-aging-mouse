library(ggplot2)
library(gridExtra)
setwd("../../analysis/snapATAC/all66")
sample_list = read.csv("../../../scripts/snapATAC/doublet_thresholds_manual.csv")
rownames(sample_list) = sample_list$Sample
outF = "all21"
load(paste(outF, ".cluster.RData",sep=""))


cluster_colors=c(
  "grey", "#E31A1C", "#FFD700", "#771122", "#777711", "#1F78B4", "#68228B", "#AAAA44",
  "#60CC52", "#771155", "#DDDD77", "#774411", "#AA7744", "#AA4455", "#117744", 
  "#000080", "#44AA77", "#AA4488", "#DDAA77", "#D9D9D9", "#BC80BD", "#FFED6F",
  "#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17",
  "#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02",
  "#A6761D", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
  "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#B15928", "#FBB4AE", "#B3CDE3",
  "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2",
  "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC",
  "#CCCCCC", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FFFF33", "#A65628",
  "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
  "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
  "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5"
)

x.sp@metaData$tissue = ""
x.sp@metaData$sample_name = ""
for(n in sample_list$Sample) {
  n = paste(n)
  x.sp@metaData[which(paste(x.sp@metaData$sample) == n),"tissue"] = paste(sample_list[n,"tissue"])
  samp_name = paste(sample_list[n,"stage"], sample_list[n,"tissue"], sample_list[n,"rep"], sep="_")
  x.sp@metaData[which(paste(x.sp@metaData$sample) == n),"sample_name"] = samp_name
}


cluster = as.character(x.sp@cluster)
cluster[which(cluster %in% c(13))]="CA1GL"
cluster[which(cluster %in% c(31))]="CA3GL"
cluster[which(cluster %in% c(12))]="PVGA"
cluster[which(cluster %in% c(15))]="ITHGL"
cluster[which(cluster %in% c(14))]="L6bGL"
cluster[which(cluster %in% c(21,22))]="ITL5GL"
cluster[which(cluster %in% c(10,18,23))]="OLG"
cluster[which(cluster %in% c(1))]="OPC"
cluster[which(cluster %in% c(8))]="MCG"
cluster[which(cluster %in% c(28,29))]="D2MSN"

cluster[which(cluster %in% c(2,5,11))]="DG"
cluster[which(cluster %in% c(4,20,30))]="AST"
cluster[which(cluster %in% c(26))]="PER"
cluster[which(cluster %in% c(24))]="VEC"
cluster[which(cluster %in% c(7))]="VLMC"
cluster[which(cluster %in% c(16))]="VPIA"
cluster[which(cluster %in% c(25))]="U1"
cluster[which(cluster %in% c(9))]="U2"
cluster[which(cluster %in% c(19))]="U3"
cluster[which(cluster %in% c(17))]="U4"


cluster = cluster[-which(cluster == 6)]
cluster = cluster[-which(cluster == 3)]
cluster = cluster[-which(cluster == 27)]


table(cluster)
x.sp = x.sp[-which(x.sp@cluster==6)]
x.sp = x.sp[-which(x.sp@cluster==3)]
x.sp = x.sp[-which(x.sp@cluster==27)]


tissue = x.sp@metaData$tissue
tissue[which(tissue %in% c("3F","4E"))] = "3F+4E"
tissue[which(tissue %in% c("5E", "6E"))] = "5E+6E"
tissue[which(tissue %in% c("7H","8H", "9G"))] = "7H+8H+9G"
tissue[which(tissue %in% c("13D","14C"))] = "13D+14C"
tissue[which(tissue %in% c("11E","11F","12E"))] ="11E+11F+12E"
tissue[which(tissue %in% c("8E","9H","8J","9J"))] = "8E+9H+8J+9J"

region =  x.sp@metaData$tissue
region[which(region %in% c("3F","4E", "3F+4E"))] = "Nucleus accumbens"
region[which(region %in% c("5E", "6E","5E+6E"))] = "Caudate Putamen"
region[which(region %in% c("7H","8H", "9G", "7H+8H+9G"))] = "Amygdala"
region[which(region %in% c("13D","14C","13D+14C"))] = "PAG/PCG"
region[which(region %in% c("11E","11F","12E","11E+11F+12E"))] ="Posterior Hippocampus"
region[which(region %in% c("8E","9H","8J","9J","8E+9H+8J+9J"))] = "Anterior Hippocampus"
region[which(region %in% c("2A+3A"))] = "Frontal Cortex"
region[which(region %in% c("12D+13B"))] = "Entorhinal Cortex"
table(region)



pdf(paste(outF, ".figure1.pdf",sep=""), height = 6, width = 18)
par(mfrow = c(1,3))
plotViz(
  obj= x.sp,
  method="umap",
  main="Cluster",
  point.color=cluster,
  point.size=0.3,
  point.shape=19,
  text.add=T,
  text.size=1,
  text.color="black",
  down.sample=10000,
  legend.add=FALSE
)

plotViz(
  obj= x.sp,
  method="umap",
  main="Brain Region",
  point.color=region,
  point.size=.3,
  point.shape=19,
  text.add=F,
  text.size=1, point.alpha = .7,
  text.color="black",
  down.sample=10000,
  legend.add=F
)


plotViz(
  obj= x.sp,
  method="umap",
  main="Age",
  point.color=x.sp@metaData$stage,
  point.size=.3,
  point.shape=19,
  text.add=F,point.alpha = .7,
  text.size=1,
  text.color="black",
  down.sample=10000,
  legend.add=F
)

par(mfrow = c(1,1))
df = data.frame(table(cluster))
colnames(df) = c("Cluster", "count")
cp<-ggplot(df, aes(x=Cluster, y=count, fill = Cluster)) +
  geom_bar(stat="identity")+scale_fill_manual(values=cluster_colors) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  text = element_text(size=20),legend.position = "none")


df = data.frame(table(region))
colnames(df) = c("Region", "count")
rp<-ggplot(df, aes(x=Region, y=count, fill = Region)) +
  geom_bar(stat="identity")+scale_fill_manual(values=cluster_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  text = element_text(size=20), legend.position = "none")


df = data.frame(table(x.sp@metaData$stage))
colnames(df) = c("Age", "count")
df$Age = factor(df$Age, levels = c("8wk","9mo","18mo"))
ap<-ggplot(df, aes(x=Age, y=count, fill = Age)) +
  geom_bar(stat="identity")+scale_fill_manual(values=cluster_colors[c(2,3,1)]) + 
  theme(legend.position = "none",  text = element_text(size=20))

print(grid.arrange(cp, rp,ap, ncol=3))

df = data.frame(table(x.sp@metaData$sample_name))
colnames(df) = c("Sample", "count")
sp<-ggplot(df, aes(x=Sample, y=count, fill = Sample)) +
  geom_bar(stat="identity")+scale_fill_manual(values=cluster_colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none")
print(sp)

par(mfrow = c(1,3))

plotViz(
  obj= x.sp,
  method="umap",
  main="brain region",
  point.color=region,
  point.size=0,
  point.shape=19,
  text.add=F,
  text.size=1,
  text.color="black",
  down.sample=10000,
  legend.add=T
)



plotViz(
  obj= x.sp,
  method="umap",
  main="Age",
  point.color=x.sp@metaData$stage,
  point.size=0,
  point.shape=19,
  text.add=F,
  text.size=1,
  text.color="black",
  down.sample=10000,
  legend.add=T
)

plotViz(
  obj= x.sp,
  method="umap",
  main="Sample Name",
  point.color=x.sp@metaData$sample_name,
  point.size=.3,
  point.shape=19,
  text.add=F,
  text.size=1,
  text.color="black",
  down.sample=10000,
  legend.add=F
)


dev.off()







#dev.off()
####





x.sp@metaData$region = x.sp@metaData$tissue

tab = data.frame(sample=x.sp@metaData$sample_name,cluster=x.sp@cluster,region = x.sp@metaData$tissue)
summ = plyr::count(tab[,1])
tab2 = plyr::count(tab)
tab2$total = summ$freq[match(tab2$sample,summ$x)]
tab2$frac = tab2$freq/tab2$total

pdf("celltype_frac_per_region.pdf",height=5,width=10)
for(i in unique(x.sp@metaData$tissue)) {
  p = ggplot(tab2[which(tab2$region==i),]) + 
    geom_col(aes(x=cluster,y=frac,fill=sample),position="dodge") +
    theme_bw() +ggtitle(i)
  print(p)
}
dev.off()




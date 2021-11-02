files = list.files("age_diff_edgeR.snap/up_bed",pattern = paste("8wk_vs_18mo"),full.names = T)
files2 = list.files("age_diff_edgeR.snap/up_bed",pattern = paste("9mo_vs_18mo"),full.names = T)
files3 = list.files("age_diff_edgeR.snap/up_bed",pattern = paste("8wk_vs_9mo"),full.names = T)
down_files = c(files,files2,files3)

files = list.files("age_diff_edgeR.snap/up_bed",pattern = paste("OLG"),full.names = T)
files1 = list.files("age_diff_edgeR.snap/down_bed",pattern = paste("OLG"),full.names = T)
olg_files = c(files,files1)

peaks = c()
for(f in files) {
  tmp = fread(f, sep = "\t")
  tpeaks = paste(tmp$V1,":",tmp$V2,"-",tmp$V3, sep = "")
  peaks = c(peaks, tpeaks)
}
length(peaks)
peaks = unique(peaks)


c = read.csv("cpm_tables/OLG_cpm.txt")
colnames(c) = gsub("X", "",colnames(c))
colnames(c)[1] = "X"
i = 2
j = 3
pdf("OLG_cluster_correlation_plots.pdf")
for (i in seq(2,48,2)) {
j = i+1
cu = c[,c(i,j)]
rownames(cu) = c$X
data = data.frame(cu)
data = data[which(rownames(data) %in% peaks),]
colors = rep("gray", nrow(data))

colnames(data) = c("X1", "X2")
g=ggplot(data %>% arrange((colors)), aes(x=X1, y=X2)) + theme(axis.text=element_text(size=14,face="bold"),axis.title=element_text(size=14,face="bold")) +
  geom_point(alpha = .5) + geom_density2d() +
  geom_abline(intercept = 0, slope = 1, colour = "green") + ggtitle(paste("CPM OLG diff peaks\n ",round(cor(data)[1,2],digits = 6), sep = "")) +
  xlab(colnames(c)[i]) + ylab(colnames(c)[j]) 
print(g)

cat(colnames(c)[c(i,j)],cor(data)[1,2],sep = "\n")
}
dev.off()

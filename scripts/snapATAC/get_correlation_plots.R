
files = list.files("age_diff_edgeR.snap/up_bed",pattern = paste("8wk_vs_18mo"),full.names = T)
files2 = list.files("age_diff_edgeR.snap/up_bed",pattern = paste("9mo_vs_18mo"),full.names = T)
files3 = list.files("age_diff_edgeR.snap/up_bed",pattern = paste("8wk_vs_9mo"),full.names = T)
down_files = c(files,files2,files3)

peaks = c()
for(f in files) {
  tmp = fread(f, sep = "\t")
  tpeaks = paste(tmp$V1,":",tmp$V2,"-",tmp$V3, sep = "")
  peaks = c(peaks, tpeaks)
}


tissue = x.sp@metaData$tissue
tissue[which(tissue %in% c("3F","4E"))] = "3F+4E"
tissue[which(tissue %in% c("5E", "6E"))] = "5E+6E"
tissue[which(tissue %in% c("7H","8H", "9G"))] = "7H+8H+9G"
tissue[which(tissue %in% c("13D","14C"))] = "13D+14C"
tissue[which(tissue %in% c("11E","11F","12E"))] ="11E+11F+12E"
tissue[which(tissue %in% c("8E","9H","8J","9J"))] = "8E+9H+8J+9J"




q = quantile(c, c(.1,.9))

smc = c[which(rowMeans(c)<as.numeric(q[2])),]
smc = c[which(rowMeans(smc)>as.numeric(q[1])),]


rlb = (cbind(rowMeans(c[,c(5,6)]),rowMeans(c[,c(51,52)])))


rlb = (cbind(c[,1],c[,2]))

library(PerformanceAnalytics)

ent = c[,grep("12D", colnames(c))]
q = quantile(rowMeans(ent), c(.4,.98))
ent = ent[-which(rowMeans(ent)>as.numeric(q[2])),]
ent = ent[-which(rowMeans(ent)<as.numeric(q[1])),]
chart.Correlation(ent[sample.int(nrow(ent),50000),])



ent = c[,grep("2A", colnames(c))]
q = quantile(rowMeans(ent), c(.4,.98))
ent = ent[-which(rowMeans(ent)>as.numeric(q[2])),]
ent = ent[-which(rowMeans(ent)<as.numeric(q[1])),]
chart.Correlation(ent[sample.int(nrow(ent),50000),])


ent = c[,unique(c(grep("8E", colnames(c)), grep("8J", colnames(c)),grep("9H", colnames(c)),grep("9J", colnames(c))))]
q = quantile(rowMeans(ent), c(.4,.98))
ent = ent[-which(rowMeans(ent)>as.numeric(q[2])),]
ent = ent[-which(rowMeans(ent)<as.numeric(q[1])),]
chart.Correlation(ent[sample.int(nrow(ent),1000),])

wk8 = mat[,grep("8wk", colnames(mat))]


wk8 = wk8[,c(grep("8E", colnames(wk8)), grep("8J", colnames(wk8)),grep("9H", colnames(wk8)),grep("9J", colnames(wk8)))]
wk8_rep1 = rowSums(wk8[,grep("rep1", colnames(wk8))])
wk8_rep2 = rowSums(wk8[,grep("rep2", colnames(wk8))])
dt8wk = data.frame(`8wk_8E+9H+8J+9J_rep1` = wk8_rep1, `8wk_8E+9H+8J+9J_rep2` = wk8_rep2 )
cdt8wk = cpm(dt8wk)
ent = cbind(c[,grep('8E+9H+8J+9J' ,colnames(c),fixed = T)],cdt8wk)
q = quantile(rowMeans(ent), c(.4,.98))
ent = ent[-which(rowMeans(ent)>as.numeric(q[2])),]
ent = ent[-which(rowMeans(ent)<as.numeric(q[1])),]
chart.Correlation(ent[sample.int(nrow(ent),50000),])




q = quantile(rowMeans(ent), c(.4,.98))
ent = ent[-which(rowMeans(ent)>as.numeric(q[2])),]
ent = ent[-which(rowMeans(ent)<as.numeric(q[1])),]
chart.Correlation(ent[peaks[which(peaks %in% rownames(ent))],])


r = cor( ent)
ggcorrplot(r, 
           hc.order = T, 
           type = "lower",
           lab = TRUE)


r = cor(c[,grep("_2A.3A_r",colnames(c))])
ggcorrplot(r, 
          # hc.order = T, 
          # type = "lower",
           lab = TRUE)
c[,c(48,49)]

data = data.frame(rlb)
rownames(data) = rownames(c)

data = data.frame(cu)
q = quantile(rowMeans(data), c(.4,.94))
data = data[-which(rowMeans(data)>as.numeric(q[2])),]
data = data[-which(rowMeans(data)<as.numeric(q[1])),]
nrow(data)
colnames(data) 
colors = rep("gray", nrow(data))
#colors[which(diff_ATAC$DLD1.HFSUMO3.vs..HFSUMO3.senp5KO.adj..p.value<0.05)] = "red"
#sam = data[which(colors == "red"),]
#show(head(sam[order(sam[,1], decreasing = T)[1:10],]))
g=ggplot(data %>% arrange((colors)), aes(x=X1, y=X2)) + theme(axis.text=element_text(size=14,face="bold"),axis.title=element_text(size=14,face="bold")) +
  geom_point(colour = colors[order(colors)]) + geom_density2d() +
  geom_abline(intercept = 0, slope = 1, colour = "green") + ggtitle(paste("CPM", sep = "")) +
  xlab(colnames(c)[1]) + ylab(colnames(c)[2])
g




g=ggplot(data %>% arrange((colors)), aes(x=X8wk_11E.11F.12E_rep1, y=X8wk_11E.11F.12E_rep2)) + theme(axis.text=element_text(size=14,face="bold"),axis.title=element_text(size=14,face="bold")) +
  geom_point(colour = colors[order(colors)]) + geom_density2d() +
  geom_abline(intercept = 0, slope = 1, colour = "green") + ggtitle(paste("CPM", sep = "")) +
  xlab(colnames(c)[1]) + ylab(colnames(c)[2])
g
nrow(data)

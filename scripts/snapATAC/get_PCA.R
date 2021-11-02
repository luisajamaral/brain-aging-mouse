


pmat = x.sp@pmat
sample = x.sp@metaData$sample_name
nsample = length(unique(sample))

  dat = list()
  for (s in unique(sample)) {
    idx1 = which(sample == s)
    dat[[s]] = colSums(pmat[idx1, ])
  }
  mat = do.call(cbind, dat)
  
  write.table(mat, file = "sample_count_mat.txt", sep = "\t", quote = F)
rm(pmat)  
head(mat)
rownames(mat) = x.sp@peak$name
library(edgeR)
c = cpm(mat)

write.table(c, file = "sample_cpm_mat.txt", sep = "\t", quote = F)

library(ggfortify)

pr <- prcomp(t(c), center = TRUE, scale = TRUE)
summary(pr)

screeplot(pr, type = "l", npcs = 15, main = "Screenplot of the first 15 PCs")
#abline(h = 1, col="red", lty=5)
#legend("topright", legend=c("Eigenvalue = 1"),
#       col=c("red"), lty=5, cex=0.6)

cumpro <- cumsum(pr$sdev^2 / sum(pr$sdev^2))
plot(cumpro[0:15], xlab = "PC #", ylab = "Amount of explained variance", main = "Cumulative variance plot")



plot(pr$x[,1],pr$x[,2], xlab="PC1 (0.1857%)", ylab = "PC2 (0.1563%)", main = "PC1 / PC2 - plot")

m = cbind(pr$x[,1],pr$x[,2])
m = as.data.frame(m)
rep = sapply(strsplit(as.character(rownames(m)), "_"), "[[", 3)
region = sapply(strsplit(as.character(rownames(m)), "_"), "[[", 2)
age = sapply(strsplit(as.character(rownames(m)), "_"), "[[", 1)

colnames(m) = c("PC1(18.57%)", "PC2(15.63%)")
ggplot(m, aes(x=`PC1(18.57%)`, y=`PC2(15.63%)`, color = region, shape = age)) +
  geom_point(size=2)
ggplot(m, aes(x=`PC1(18.57%)`, y=`PC2(15.63%)`, color = age)) +
  geom_point(size=2)
ggplot(m, aes(x=`PC1(18.57%)`, y=`PC2(15.63%)`, color = rep)) +
  geom_point(size=2)
ggplot(m, aes(x=`PC1(18.57%)`, y=`PC2(15.63%)`, color = region)) +
  geom_point(size=2)

dev.off()











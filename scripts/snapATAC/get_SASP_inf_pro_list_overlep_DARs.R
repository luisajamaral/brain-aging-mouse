

# Basic function to convert human to mouse gene names
convertHumanGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}


SASP = read.csv("SASP.csv",header = F)
SASP = as.character(SASP$V1)
SASP_genes <- convertHumanGeneList(SASP)

pro = read.csv("proliferation.csv",header = F)
pro = as.character(pro$V1)
pro_genes <- convertHumanGeneList(pro)

inf = read.csv("inflammatory.csv",header = F)
inf = as.character(inf$V1)
inf_genes <- convertHumanGeneList(inf)


genes = read.table("../../../annotations/mm10.gencode.vM16.gene.bed")
# looking for all peaks in and around the gene (not just promoters)
genes.gr = GRanges(genes[,1],
                   IRanges(genes[,2]-2000, genes[,3]+2000), name=genes[,4]
)
s_granges = genes.gr[which(genes.gr$name %in% c(SASP_genes,pro_genes,inf_genes))]
st = data.frame(cbind(c(SASP_genes,pro_genes,inf_genes), c(rep("SASP", length(SASP_genes)), rep("proliferation",length(pro_genes)), rep("inflammation", length(inf_genes)))))
st$type = ""

for(i in unique(st$X1)) {
  st$type[which(st$X1==i)] = paste(st[which(st$X1== i),"X2"], collapse = ",")
}
st = st[-which(duplicated(st$X1)),]

rownames(st) = st$X1
st = as.data.frame(st)
st = st[,-c(2)]
st


rename_file = function(file) {
  file = gsub("13","CA1GL", file)
  file = gsub("31","CA3GL", file)
  file = gsub("12","PVGA", file)
  file = gsub("15","ITHGL", file)
  file = gsub("14","L6bGL", file)
  file = gsub("21","ITL5GL", file)
  file = gsub("22","ITL5GL2", file)
  file = gsub("26","PER", file)
  file = gsub("24","VEC", file)
  file = gsub("^7","VLMC", file)
  file = gsub("16","VPIA", file)
  return(file)
}

files = list.files("age_diff_edgeR.snap/", "9mo_vs_18mo")
for (file in files) {
  d = read.csv(paste("age_diff_edgeR.snap/",file,sep = ""))
  cat(nrow(d[which(d$fdr<0.1),]))
  file = rename_file(file)
  coors = d$X
  chr = sapply(strsplit(as.character(coors), ":"), "[[", 1)
  coors = sapply(strsplit(as.character(coors), ":"), "[[", 2)
  start = sapply(strsplit(as.character(coors), "-"), "[[", 1)
  end = sapply(strsplit(as.character(coors), "-"), "[[", 2)
  d_granges = GRanges(chr,
                      IRanges(as.numeric(start), as.numeric(end)), fdr=d$fdr, logFC = -d$logFC
  )
  
  
  f  = findOverlaps(d_granges,s_granges)
  result = cbind(paste(decode(d_granges@seqnames)[f@from],d_granges[f@from]@ranges, sep = ":"),
                 paste(d_granges[f@from]$logFC),paste(d_granges[f@from]$fdr),
                 paste(s_granges[f@to]$name),paste(decode(s_granges@seqnames)[f@to],
                                                   s_granges[f@to]@ranges,sep = ":"))
  result = data.frame(result)
  result$type = st[paste(result$X4),2]
  result[1:15,]
  result = as.data.frame(result)
  result = result[-which(duplicated(result$X4)),]
  result = result[which(as.numeric(as.character(result$X3)) < 0.05),]
  colnames(result) = c("peak","logFC","FDR","gene","gene_region","type")
  rownames(result) = result$gene
  
  result = result[,-c(4)]
  show(head(result))
  if(nrow(result)>=1) {
  write.table(result, paste("SASP_pro_inf_overlap/",gsub(".edger.txt","",file),"_overlap.txt",sep = ""), sep = "\t", quote = F)
  }
}

files = list.files("age_diff_edgeR.snap/", "8wk_vs_18mo")
sink("SASP_pro_inf_overlap/8wk_vs_18mo/num_DAR_total.txt")
for (file in files) {
  d = fread(paste("age_diff_edgeR.snap/",file,sep = ""))
  file = rename_file(file)
  
  cat(file, nrow(d[which(d$fdr<0.05),]), "\n")
}
sink()

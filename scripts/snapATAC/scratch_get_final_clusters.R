
# re annotate final 

load("../../analysis/snapATAC/all66/all.cluster.meta.CEMBAann.glial.RData")
glial = x.sp.sub

glial_meta = glial@metaData
rownames(glial_meta) = paste(glial_meta$sample,glial_meta$barcode, sep = "_")

load("../../analysis/snapATAC/all66/all.cluster.meta.CEMBAann.AST.RData")
x.sp.sub@metaData$cluster_names_final = "AST"
x.sp.sub@metaData$cluster_names_final[which(x.sp.sub@cluster == 3)]= "RGDG/NIPC"
x.sp.sub@metaData$cluster_names_final[which(x.sp.sub@cluster == 16)]= "RGL"


ast_meta = x.sp.sub@metaData
rownames(ast_meta) = paste(ast_meta$sample,ast_meta$barcode, sep = "_")
r = rownames(ast_meta[which(ast_meta$cluster_names_final=="RGL"),])
glial_meta[r,"cluster_names_final"] = "RGL"
r = rownames(ast_meta[which(ast_meta$cluster_names_final=="RGDG/NIPC"),])
glial_meta[r,"cluster_names_final"] = "RGDG/NIPC"



load("../../analysis/snapATAC/all66/all.cluster.meta.CEMBAann.ExN.RData")
ExN_meta =x.sp.sub@metaData
rownames(ExN_meta) = paste(ExN_meta$sample,ExN_meta$barcode, sep = "_")

load("../../analysis/snapATAC/all66/all.cluster.meta.CEMBAann.InhN.RData")
InhN_meta =x.sp.sub@metaData
rownames(InhN_meta) = paste(InhN_meta$sample,InhN_meta$barcode, sep = "_")

all_meta = rbind(glial_meta, ExN_meta, InhN_meta)
load("../../analysis/snapATAC/all66/all.cluster.meta.final.RData")

m = x.sp@metaData
rownames(m) = paste(m$sample,m$barcode, sep = "_")

x.sp@metaData$final_clusters = all_meta[paste(rownames(m)),"cluster_names_final"]
x.sp@metaData$final_clusters[which(x.sp@metaData$final_clusters=="Amy_2")] = "Amy_InhN"
x.sp@metaData$final_clusters[which(x.sp@metaData$final_clusters=="Amy_1")] = "Amy_InhN"
x.sp@metaData$final_clusters[which(x.sp@metaData$final_clusters=="Amy_neur")] = "Amy_ExN"

x.sp@metaData$final_clusters[which(x.sp@metaData$final_clusters=="PAG_n")] = "PAG/PCG_N"

pdf("../../analysis/snapATAC/all66/final_cluster.UMAP.pdf")
plotViz(
  obj= x.sp,
  method="umap",
  main="Final annotation",
  point.size=0.3,
  point.shape=19,
  point.color=x.sp@metaData$final_clusters,
  text.add=T,
  text.size=.5,
  text.color="black",
  down.sample=10000,
  legend.add=F
)
dev.off()


save(x.sp,file= "../../analysis/snapATAC/all66/all.cluster.meta.final.RData")


x.sp@metaData$final_clusters[which(x.sp@metaData$final_clusters %in% c("D1MSN","D2MSN"))] = "MSN"

load("../../analysis/snapATAC/all66/all.cluster.meta.CEMBAann.AST.RData")
x.sp.sub@metaData$cluster_names_final = "AST"
x.sp.sub@metaData$cluster_names_final[which(x.sp.sub@cluster == 3)]= "RGDG/NIPC"
x.sp.sub@metaData$cluster_names_final[which(x.sp.sub@cluster == 16)]= "RGL"

x.sp@metaData$final_clusters[which(x.sp@metaData$cluster_names=="AST")] = "AST"
rownames(x.sp@metaData) = paste(x.sp@metaData$sample,x.sp@metaData$barcode, sep = "_")


ast_meta = x.sp.sub@metaData
rownames(ast_meta) = paste(ast_meta$sample,ast_meta$barcode, sep = "_")
r = rownames(ast_meta[which(ast_meta$cluster_names_final=="RGL"),])
x.sp@metaData[r,"final_clusters"] = "RGL"
r = rownames(ast_meta[which(ast_meta$cluster_names_final=="RGDG/NIPC"),])
x.sp@metaData[r,"final_clusters"] = "RGDG/NIPC"



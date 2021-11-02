library("leiden")
load("../../analysis/snapATAC/all66/all.cluster.meta.CEMBAann.glial.RData")
tic("runCluster_R-igraph")
x.sp.sub = runCluster(obj=x.sp.sub, 
                      tmp.folder=tempdir(), 
                      louvain.lib = "leiden",
                      seed.use=10,resolution = 1.5
);
toc()


plotViz(
  obj= x.sp.sub,
  method="umap",
  main="Cluster",
  point.color=x.sp.sub@cluster,
  point.size=0.2,
  point.shape=19,
  text.add=T,
  text.size=1,
  text.color="black",
  down.sample=10000,
  legend.add=FALSE
)

ca = paste( x.sp.sub@metaData$CEMBA_major,x.sp.sub@cluster)
table(ca)
# x.sp.sub@metaData$cluster_names_final = ""#x.sp.sub@metaData$cluster_names
# x.sp.sub@metaData$cluster_names_final[which(x.sp.sub@cluster %in% c(1,2,3))] = "GRC"
# x.sp.sub@metaData$cluster_names_final[which(x.sp.sub@cluster %in% c(6))] = "CA3GL"
# x.sp.sub@metaData$cluster_names_final[which(x.sp.sub@cluster %in% c(5,7))] = "CA1GL"
# x.sp.sub@metaData$cluster_names_final[which(x.sp.sub@cluster %in% c(14))] = "NPGL"
# x.sp.sub@metaData$cluster_names_final[which(x.sp.sub@cluster %in% c(12,16))] = "ITHGL"
# x.sp.sub@metaData$cluster_names_final[which(x.sp.sub@cluster %in% c(13))] = "CTGL"
# x.sp.sub@metaData$cluster_names_final[which(x.sp.sub@cluster %in% c(11))] = "L6bGL"
# x.sp.sub@metaData$cluster_names_final[which(x.sp.sub@cluster %in% c(4))] = "ITL6GL"
# x.sp.sub@metaData$cluster_names_final[which(x.sp.sub@cluster %in% c(8))] = "ITL23GL"
# x.sp.sub@metaData$cluster_names_final[which(x.sp.sub@cluster %in% c(10))] = "Ent_neur"
# x.sp.sub@metaData$cluster_names_final[which(x.sp.sub@cluster %in% c(9))] = "Amy_neur"
# x.sp.sub@metaData$cluster_names_final[which(x.sp.sub@cluster %in% c(15))] = "CA3GL"



x.sp.sub@metaData$cluster_names_final = ""#x.sp.sub@metaData$cluster_names
x.sp.sub@metaData$cluster_names_final[which(x.sp.sub@cluster %in% c(1))] = "PAG_n"
x.sp.sub@metaData$cluster_names_final[which(x.sp.sub@cluster %in% c(4))] = "Amy_1"
x.sp.sub@metaData$cluster_names_final[which(x.sp.sub@cluster %in% c(5))] = "Amy_2"
x.sp.sub@metaData$cluster_names_final[which(x.sp.sub@cluster %in% c(8))] = "SSTGA"
x.sp.sub@metaData$cluster_names_final[which(x.sp.sub@cluster %in% c(7))] = "PVGA"
x.sp.sub@metaData$cluster_names_final[which(x.sp.sub@cluster %in% c(10))] = "VIPGA"
x.sp.sub@metaData$cluster_names_final[which(x.sp.sub@cluster %in% c(11))] = "LAMGA"
x.sp.sub@metaData$cluster_names_final[which(x.sp.sub@cluster %in% c(3))] = "D1MSN"
x.sp.sub@metaData$cluster_names_final[which(x.sp.sub@cluster %in% c(2))] = "D2MSN"
x.sp.sub@metaData$cluster_names_final[which(x.sp.sub@cluster %in% c(6))] = "STRGA2"
x.sp.sub@metaData$cluster_names_final[which(x.sp.sub@cluster %in% c(9))] = "STRGA1"




x.sp.sub@metaData$cluster_names_final = x.sp.sub@metaData$cluster_names
x.sp.sub@metaData$cluster_names_final[which(x.sp.sub@cluster %in% c(14))] = "IOL"


x.sp.sub@metaData$cluster_names_final = "AST"
x.sp.sub@metaData$cluster_names_final[which(x.sp.sub@cluster %in% c(10))] = "RGL"



pdf("Ast_final_clusters.pdf")
plotViz(
  obj= x.sp.sub,
  method="umap",
  main="Cluster",
  point.color=x.sp.sub@metaData$cluster_names_final,
  point.size=0.2,
  point.shape=19,
  text.add=T,
  text.size=1,
  text.color="black",
  down.sample=10000,
  legend.add=F
)
dev.off()

save(x.sp.sub,file = "../../analysis/snapATAC/all66/all.cluster.meta.CEMBAann.Ast.RData")


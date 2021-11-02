
t = fread("/projects/ps-renlab/yangli/projects/CEMBA/01.joint_dat/rs1cemba/rs1cemba.full.cellMeta.taxonomyLable.tsv")
b = paste0(t$Barcode,sapply(strsplit(as.character(t$Sample), "_"), "[[", 2))
ky = unique(t$Sample)
ky = ky[order(sapply(strsplit(as.character(ky), "_"), "[[", 2))]
ky_orig = ky
ky[seq(1,length(ky),2)] = paste(ky[seq(1,length(ky),2)],"rep1",sep = "_")
ky[seq(2,length(ky),2)] = paste(ky[seq(2,length(ky),2)],"rep2",sep = "_")
ky = paste(sapply(strsplit(as.character(ky), "_"), "[[", 2), sapply(strsplit(as.character(ky), "_"), "[[", 3),sep = "_")
dt = data.frame(cbind(ky_orig,ky))
rownames(dt) = dt$ky_orig
t$cell.id = paste(dt[t$Sample,2],t$Barcode, sep = "_")
t$samp = paste(dt[t$Sample,2])


load("../../analysis/snapATAC/all66/sub_ExN.cluster.RData")
x.sp.sub@metaData$cell.id = x.sp.sub@metaData$sample_name
x.sp.sub@metaData$cell.id = paste(sapply(strsplit(as.character(x.sp.sub@metaData$cell.id), "_"), "[[", 2), sapply(strsplit(as.character(x.sp.sub@metaData$cell.id), "_"), "[[", 3),sep = "_")
x.samp = x.sp.sub@metaData$cell.id

x.sp.sub@metaData$cell.id = paste(x.sp.sub@metaData$cell.id,x.sp.sub@metaData$barcode,sep = "_")
x.sp.sub@metaData$cell.id[grep("9mo",x.sp.sub@metaData$sample_name)] = ""
x.sp.sub@metaData$cell.id[grep("18mo",x.sp.sub@metaData$sample_name)] = ""

t = t[which(t$samp %in% x.samp),]


x.sp.sub@metaData$CEMBA_major = ""
x.sp.sub@metaData$CEMBA_subtype = ""
x.sp.sub@metaData$class = ""

t = as.data.frame(t)
rownames(t) = t$cell.id
o = t[paste( x.sp.sub@metaData$cell.id),]


x.sp.sub@metaData$CEMBA_major = o$MajorType
x.sp.sub@metaData$CEMBA_subtype = o$SubType

o$MajorType = x.sp.sub@metaData$CEMBA_major
o$MajorType[which(o$MajorType %in% names(which(table(o$MajorType)<15)))] = NA

o$SubType = x.sp.sub@metaData$CEMBA_subtype
o$SubType[which(o$SubType %in% names(which(table(o$SubType)<15)))] = NA



pdf("../../analysis/snapATAC/all66/CEMBA_color_ann_UMAP_AST.pdf")
plotViz(
  obj= x.sp.sub,
  method="umap",
  main="CEMBA annotation",
  point.size=0.3,
  point.shape=19,
  point.color=x.sp.sub@cluster,
  text.add=T,
  text.size=.8,
  text.color="black",
  #down.sample=10000,
  legend.add=F
)
plotViz(
  obj= x.sp.sub[-which(is.na(x.sp.sub@metaData$CEMBA_major)),],
  method="umap",
  main="CEMBA annotation",
  point.size=0.3,
  point.shape=19,
  point.color=o$MajorType[-which(is.na(x.sp.sub@metaData$CEMBA_major))],
  text.add=T,
  text.size=.8,
  text.color="black",
  #down.sample=10000,
  legend.add=F
)
plotViz(
  obj= x.sp.sub[-which(is.na(x.sp.sub@metaData$CEMBA_major)),],
  method="umap",
  main="CEMBA annotation",
  point.size=0.3,
  point.shape=19,
  point.color=o$MajorType[-which(is.na(x.sp.sub@metaData$CEMBA_major))],
  text.add=FALSE,
  text.size=1.5,
  text.color="black",
  #down.sample=10000,
  legend.add=F
)
plotViz(
  obj= x.sp.sub[-which(is.na(x.sp.sub@metaData$CEMBA_major)),],
  method="umap",
  main="CEMBA annotation",
  point.size=0,
  point.shape=19,
  point.color=o$MajorType[-which(is.na(x.sp.sub@metaData$CEMBA_major))],
  text.add=FALSE,
  text.size=1.5,
  text.color="black",
  #down.sample=10000,
  legend.add=T
)
plotViz(
  obj= x.sp.sub[-which(is.na(x.sp.sub@metaData$CEMBA_subtype)),],
  method="umap",
  main="CEMBA annotation",
  point.size=0.3,
  point.shape=19,
  point.color=o$SubType[-which(is.na(x.sp.sub@metaData$CEMBA_subtype))],
  text.add=F,
  text.size=.7,
  text.color="black",
  #down.sample=10000,
  legend.add=F
)
plotViz(
  obj= x.sp.sub[-which(is.na(x.sp.sub@metaData$CEMBA_subtype)),],
  method="umap",
  main="CEMBA annotation",
  point.size=0,
  point.shape=19,
  point.color=o$SubType[-which(is.na(x.sp.sub@metaData$CEMBA_subtype))],
  text.add=F,
  text.size=.7,
  text.color="black",
  #down.sample=10000,
  legend.add=T
)
plotViz(
  obj= x.sp.sub[-which(is.na(x.sp.sub@metaData$CEMBA_subtype)),],
  method="umap",
  main="CEMBA annotation",
  point.size=0.3,
  point.shape=19,
  point.color=o$SubType[-which(is.na(x.sp.sub@metaData$CEMBA_subtype))],
  text.add=T,
  text.size=.7,
  text.color="black",
 # down.sample=10000,
  legend.add=F
)
dev.off()

pdf("../../analysis/snapATAC/all66/CEMBA_color_ann_UMAP_split_AST.pdf")
u = unique(o$SubType)
u = u[-which(is.na(u))]
par(mfrow = c(3,2))
for (ty in u) {
  plotViz(
    obj= x.sp.sub[which(x.sp.sub@metaData$CEMBA_subtype == ty),],
    method="umap",
    main=ty,
    point.size=0.3,
    point.shape=19,
    point.color=o$MajorType[which(x.sp.sub@metaData$CEMBA_subtype == ty)],
    text.add=F,
    text.size=1,
    text.color="black",
    #down.sample=10000,
    legend.add=F
  )
}
dev.off()

save(x.sp.sub,file = "../../analysis/snapATAC/all66/all.cluster.meta.CEMBAann.D2MSN.RData")





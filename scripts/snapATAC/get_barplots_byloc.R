o = fread("../OLG_cpm.txt", sep = ",")
rownames(o) = o$V1
o = o[,-1]
o = data.frame(o)
o["chrY:90741087-90744992",]
region = sapply(strsplit(as.character(colnames(o)), "_"), "[[", 2)
rep = sapply(strsplit(as.character(colnames(o)), "_"), "[[", 3)
stage = sapply(strsplit(as.character(colnames(o)), "_"), "[[", 1)

loc = "chrY:90741087-90744992"
loc = "chr13:21809701-21810979"
loc = "chrY:90783703-90785266"

for (f in list.files(".","cpm")) {
  o = fread(f, sep = ",")
  ct = gsub("_cpm.txt", "", f)
  o = data.frame(o)
  rownames(o) = o$V1
  o = o[,-1]
  o["chrY:90741087-90744992",]
  region = sapply(strsplit(as.character(colnames(o)), "_"), "[[", 2)
  rep = sapply(strsplit(as.character(colnames(o)), "_"), "[[", 3)
  stage = sapply(strsplit(as.character(colnames(o)), "_"), "[[", 1)
  
  pdf(paste(ct,"barplots.pdf", sep = "_"))
  for(loc in c("chrY:90741087-90744992","chr13:21809701-21810979","chrY:90783703-90785266")) {
    d = melt(o[loc,])
    d$region = region
    d$rep = rep
    d$stage = stage
    d$stage = gsub("X","",d$stage)
    d$variable = gsub("X","",d$variable)
    d$stage = factor(d$stage,levels = c("8wk", "9mo","18mo"), ordered = T)
    
    d$stage_rep = paste(d$stage,"_",d$rep,sep = "")
    d$stage_rep = factor(d$stage_rep,levels = c("8wk_rep1","8wk_rep2" ,"9mo_rep1", "9mo_rep2","18mo_rep1", "18mo_rep2"), ordered = T)
    
   print( ggplot(d) + 
      geom_col(aes(factor(stage_rep),value,fill=factor(stage)),position="dodge")+
      xlab("Sample") + ylab("CPM")+
      theme_bw() + theme(axis.text.x = element_text(angle = 90)) + ggtitle(label = loc)+ 
      facet_wrap(~region))
    
    # require(gridExtra)
    # g = list()
    # for (r in unique(region)) {
    #   c = d[grep(r,d$variable),]
    #   p = ggplot(c) + 
    #     geom_col(aes(factor(stage_rep),value,fill=factor(stage)),position="dodge")+
    #     xlab("Sample") + ylab("CPM")+
    #     theme_bw() + theme(axis.text.x = element_text(angle = 90)) + ggtitle(label = r,subtitle = loc) 
    #   g[[r]] = p
    # }
    # print(grid.arrange(g[[1]], g[[2]],g[[3]],g[[4]],g[[5]], g[[6]],g[[7]],g[[8]], ncol = 3))
  }
  dev.off()
}

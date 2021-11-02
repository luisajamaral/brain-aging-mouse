

head(mat)
g = glm.nb(mat)
power.nb.test(mu0 = 0.8, RR = 1/.8, theta = .4, duration = 1, power = 0.8, approach = 1)


power.prop.test(p1 = .1,p2 = .2, sig.level = .05,power=.8)
effect_size = proportion_effectsize(0.7,0.8)
sample_size = zt_ind_solve_power(effect_size=effect_size, nobs1=None, alpha=0.05, power=0.8, ratio=1)


pdf("power_analysis.pdf",height = 4, width = 4)
for(frac in c(.05,.1,.2,.3,.4,.5)){
mat = matrix(0,nrow=10,ncol=10)
for (i in seq(.1,.9, .1)){
  for (j in seq(.1,.9, .1)){
    p1 = j
    p2 = min(j+frac*j, 1)
    cat(p1,p2,"\n")
    p=power.prop.test(p1 = p1,p2 =p2, sig.level = .05,power=i)
    mat[i*10,j*10] = p$n
  }
}
mat
mat = mat[c(1:9),c(1:9)]
rownames(mat) = seq(.1,.9, .1)
colnames(mat) = seq(.1,.9, .1)
print(pheatmap(mat,cluster_rows = F, cluster_cols = F,display_numbers = T,number_format = "%.1f",fontsize_number = 7,show_rownames = T,show_colnames = T,number_color = "black",
         main=paste(frac*100,"% change", sep = "")))
}
dev.off()

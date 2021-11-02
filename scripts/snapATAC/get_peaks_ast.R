library(SnapATAC)
setwd("~/projects/brain_aging_mouse/analysis/snapATAC/all66")
load(paste("all21.cluster.RData",sep=""))
runMACS(
          obj=x.sp[which(x.sp@cluster %in% c(4,20,30)),], 
          output.prefix="all21.Ast",
          path.to.snaptools="/projects/ps-renlab/lamaral/software/miniconda3/envs/py27/bin/snaptools",
          path.to.macs="/projects/ps-renlab/lamaral/software/miniconda3/envs/py27/bin/macs2",
          gsize="mm", 
          num.cores=5,
          buffer.size=500, 
          macs.options="--nomodel --shift 37 --ext 73 --qval 1e-2 -B --SPMR --call-summits",
          tmp.folder=tempdir()
          )

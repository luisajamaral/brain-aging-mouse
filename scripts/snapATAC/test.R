library("SnapATAC")
library("reticulate")
#load("../../analysis/by_sample.snapATAC/8E-9H-8J-9J_09_rep1/snapFiles/8E-9H-8J-9J_09_rep1.Frag1000.TSS10.AllCells.seed1.dimPC20.K20.res0.7.cluster.RData")
#x.sp = addPmatToSnap(x.sp, 6)
#mat.use = x.sp@pmat
#setSessionTimeLimit(cpu = Inf, elapsed = Inf)
os <- import("os", convert = FALSE)
scipy.io <- import("scipy.io" ,convert = FALSE)
import_from_path("six", path
                 = "/projects/ps-renlab/lamaral/software/miniconda3/envs/py3.6/lib/python3.6/site-packages/")
np <- import("numpy", convert = FALSE)
scr <- import("scrublet", convert = FALSE)
#mat.use = r_to_py(mat.use)



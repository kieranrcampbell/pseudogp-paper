
## Analysis of p-values

library(readr)
library(ggplot2)
library(cowplot)
library(MCMCglmm)
library(coda)
library(matrixStats)
library(reshape2)
library(scater)
library(embeddr)
library(dplyr)
library(rhdf5)



## get the paths right
#base_dir <- "~/mount/"
base_dir <- "/net/isi-scratch/kieran/"
source(file.path(base_dir, "pseudogp-paper/analysis/diffexpr/common.R"))


## Edit 'monocle' here ------
outputfile <- file.path(base_dir, "pseudogp-paper/analysis/diffexpr/waterfall/all_plots.pdf")
h5_diffexpr <- file.path(base_dir, "pseudogp-paper/data/waterfall_diffexpr.h5")
pstfile <- file.path(base_dir, "pseudogp-paper/data/waterfall_stan_traces.h5")
fdrfile <- file.path(base_dir, "pseudogp-paper/analysis/diffexpr/waterfall_fdr.txt")

load(file.path(base_dir, "pseudogp-paper/data/sce_waterfall.Rdata"))

##------- Waterfall ONLY
exprs(sce) <- log2(exprs(sce) + 1)
sce@lowerDetectionLimit <- 0.1
## end


sigList <- pvalsFromHDF5(h5_diffexpr)

pst <- h5read(pstfile, "pst")

## need sce & pst
generatePlots(sigList, sce, pst, outputfile, fdrfile)



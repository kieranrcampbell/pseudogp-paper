
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
base_dir <- "~/mount/"
source(file.path(base_dir, "pseudogp-paper/analysis/diffexpr/common.R"))


## Edit 'monocle' here ------
outputfile <- file.path(base_dir, "pseudogp-paper/analysis/diffexpr/monocle/all_plots.pdf")
h5_diffexpr <- file.path(base_dir, "pseudogp-paper/data/monocle_diffexpr.h5")
pstfile <- file.path(base_dir, "pseudogp-paper/data/monocle_stan_traces.h5")
fdrfile <- file.path(base_dir, "pseudogp-paper/analysis/diffexpr/monocle_fdr.txt")

##------- monocle ONLY
sce_file <- file.path(base_dir, "pseudogp-paper/data/sce_monocle.Rdata")
load(sce_file)
sce <- sce_23
## end

sigList <- pvalsFromHDF5(h5_diffexpr)

pst <- h5read(pstfile, "pst")

## need sce & pst
generatePlots(sigList, sce, pst, outputfile, fdrfile)



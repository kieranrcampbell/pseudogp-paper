
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

outputfile <- file.path(base_dir, "pseudogp-paper/analysis/diffexpr/ear/all_plots.pdf")
h5_diffexpr <- file.path(base_dir, "pseudogp-paper/data/ear_diffexpr.h5")
pstfile <- file.path(base_dir, "pseudogp-paper/data/ear_stan_traces.h5")
fdrfile <- file.path(base_dir, "pseudogp-paper/analysis/diffexpr/ear_fdr.txt")

load(file.path(base_dir, "pseudogp-paper/data/sce_ear.Rdata"))
sce <- sct


sigList <- pvalsFromHDF5(h5_diffexpr)

pst <- h5read(pstfile, "pst")

## need sce & pst
generatePlots(sigList, sce, pst, outputfile, fdrfile)



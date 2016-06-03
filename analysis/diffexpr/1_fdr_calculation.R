
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

  
source("analysis/diffexpr/common.R")

args <- commandArgs(trailingOnly = TRUE)
pstfile <- args[1]
h5_diffexpr <- args[2]
sce_file <- args[3]
outputfile <- args[4]
fdrfile <- args[5]


load(sce_file)

sigList <- pvalsFromHDF5(h5_diffexpr)

pst <- h5read(pstfile, "pst")

## need sce & pst
generatePlots(sigList, sce, pst, outputfile, fdrfile)

#' Performs differential gene testing for a given cell in a given representation

library(scater)
library(embeddr)
library(MCMCglmm)
library(coda)
library(rstan)
library(readr)

args <- commandArgs(trailingOnly = TRUE)

stopifnot(length(args) == 1)

i <- as.numeric(args[1]) # which posterior trace are we performing inference on?

load("data/sce_trapnell.Rdata")

load("data/resamples/gplvm_fit_all.Rdata")
#load("data/resamples/pca_resamples.Rdata")

n_cells_exprs <- rowSums(exprs(sce) > sce@lowerDetectionLimit)
genes_to_use <- n_cells_exprs > (0.1 * ncol(sce)) # select genes expressed in at least 10% of cells
sce <- sce[genes_to_use,]

extr <- extract(fit, "t")$t
sce$pseudotime <- extr[i,]
pvals <- pseudotime_test(sce, n_cores = 1)

csv_file <- paste0("data/resamples/all_cells_diffexpr/pvals_", i, ".csv")
write_csv(pvals, csv_file)

library(scater)
library(embeddr)
library(MCMCglmm)
library(coda)
library(rstan)
library(readr)

args <- commandArgs(trailingOnly = TRUE)

i <- as.numeric(args[1]) # which resample are we performing inference on?

load("data/resamples/sce_trapnell_resamples.Rdata")

load(paste0("data/resamples/gplvm_fits/fit_", i, ".Rdata"))
load("data/resamples/pca_resamples.Rdata")

which_cells <- PCA_reps[[i]]$which_cells

sce <- sce[, which_cells]
sce$pseudotime <- posterior.mode(mcmc(extract(fit, "t")$t))
pvals <- pseudotime_test(sce, n_cores = 1)

csv_file <- paste0("data/resamples/diffexpr/pvals_",i,".csv")
write_csv(pvals, csv_file)
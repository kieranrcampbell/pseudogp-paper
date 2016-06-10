
# Differential expression comparison for uncertainty in pseudotime
# kieran.campbell@sjc.ox.ac.uk

library(scater)
library(embeddr)
library(readr)
library(MCMCglmm)
library(coda)
library(rhdf5)

args <- commandArgs(trailingOnly = TRUE)
study <- args[1]
csv_file <- args[2]
stopifnot(study %in% c("trapnell", "burns", "shin"))

tracefile <- paste0("data/", study, "_pseudotime_traces.h5")
scefile <- paste0("data/sce_", study, ".Rdata")


# Load pseudotime assignments ------
if(!require(rhdf5)) stop("Need some hdf5 love")


pst <- h5read(tracefile, "pst")

load(sce_file)
sce@lowerDetectionLimit <- 0.1


## subset to genes expressed in at lesat 10% of cells
n_cells_exprs <- rowSums(exprs(sce) > sce@lowerDetectionLimit)
genes_to_use <- n_cells_exprs > (0.1 * ncol(sce)) # select genes expressed in at least 10% of cells
sce <- sce[genes_to_use,]

# Compute map estimate
tmap <- posterior.mode(mcmc(pst))

# What range of pst are we looking at?
sce$pseudotime <- tmap


ss_test <- pseudotime_test(sce, n_cores = 1)
write_csv(ss_test, csv_file)



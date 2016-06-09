
# Differential expression comparison for uncertainty in pseudotime
# kieran.campbell@sjc.ox.ac.uk

library(devtools)
library(scater)
library(embeddr)
library(switchde)
library(readr)

source("analysis/diffexpr/common.R")


tracefile <- csv_file <- sce_file <- NULL

args <- commandArgs(trailingOnly = TRUE)

if(length(args) > 0) {
  tracefile <- args[1]
  sce_file <- args[2]
  trace <- as.numeric(args[3])
  csv_file <- args[4]
} else {
  stop("Provide start and end samples")
}

# Load pseudotime assignments ------
if(!require(rhdf5)) stop("Need some hdf5 love")

print(paste("Performing differential expression analysis for trace", trace, "with csv file", csv_file))

pst <- h5read(tracefile, "pst")

load(sce_file)
sce@lowerDetectionLimit <- 0.1

## sanity checking
stopifnot(start > 0 && end > 0 && start < nrow(pst) && end < nrow(pst) && start < end)


## subset to genes expressed in at lesat 10% of cells
n_cells_exprs <- rowSums(exprs(sce) > sce@lowerDetectionLimit)
genes_to_use <- n_cells_exprs > (0.1 * ncol(sce)) # select genes expressed in at least 10% of cells
sce <- sce[genes_to_use,]

# What range of pst are we looking at?
sce$pseudotime <- pst[trace,]


ss_test <- pseudotime_test(sce, n_cores = 1)
write_csv(ss_test, output_csv)



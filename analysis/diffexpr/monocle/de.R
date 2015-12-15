
# Differential expression comparison for uncertainty in pseudotime
# kieran.campbell@sjc.ox.ac.uk

library(devtools)
library(scater)
library(embeddr)
library(switchde)

base_dir <- "/net/isi-scratch/kieran/" # needs for cluster running
source(file.path(base_dir, "pseudogp-paper/analysis/diffexpr/common.R"))

start <- NULL
end <- NULL

args <- commandArgs(trailingOnly = TRUE)

if(length(args) > 0) {
  start <- as.numeric(args[1])
  end <- as.numeric(args[2])
} else {
  stop("Provide start and end samples")
}

# Load pseudotime assignments ------
if(!require(rhdf5)) stop("Need some hdf5 love")

##--------------- Edit here
tracefile <- paste0(base_dir, "pseudogp-paper/data/monocle_stan_traces.h5")
h5outputfile <- paste0(base_dir, "pseudogp-paper/data/monocle_diffexpr.h5")
##--------------- End edit


pst <- h5read(tracefile, "pst")

##------- MONOCLE ONLY
sce_file <- file.path(base_dir, "pseudogp-paper/data/sce_monocle.Rdata")
load(sce_file)
sce <- sce_23
## end


## sanity checking
stopifnot(start > 0 && end > 0 && start < nrow(pst) && end < nrow(pst) && start < end)

doDiffExprAnalysis(sce, pst, start, end, h5outputfile)
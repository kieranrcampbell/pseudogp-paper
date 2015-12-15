
# Differential expression comparison for uncertainty in pseudotime
# kieran.campbell@sjc.ox.ac.uk

library(devtools)
library(scater)
library(embeddr)
library(switchde)
library(rhdf5)

base_dir <- "/net/isi-scratch/kieran/"
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


##--------------- Edit here
tracefile <- paste0(base_dir, "pseudogp-paper/data/waterfall_stan_traces.h5")
scefile <- paste0(base_dir, "pseudogp-paper/data/sce_waterfall.Rdata")
h5outputfile <- paste0(base_dir, "pseudogp-paper/data/waterfall_diffexpr.h5")
##--------------- End edit


load(scefile) ## loads in sce
pst <- h5read(tracefile, "pst")

##------- Waterfall ONLY
exprs(sce) <- log2(exprs(sce) + 1)
sce@lowerDetectionLimit <- 0.1
## end


## sanity checking
stopifnot(start > 0 && end > 0 && start < nrow(pst) && end < nrow(pst) && start < end)

doDiffExprAnalysis(sce, pst, start, end, h5outputfile)






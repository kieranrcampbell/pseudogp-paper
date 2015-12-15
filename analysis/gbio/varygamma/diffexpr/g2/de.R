
# Differential expression comparison for uncertainty in pseudotime
# kieran.campbell@sjc.ox.ac.uk

library(devtools)
library(scater)
library(embeddr)


start <- NULL
end <- NULL

args <- commandArgs(trailingOnly = TRUE)

if(length(args) > 0) {
  start <- as.numeric(args[1])
  end <- as.numeric(args[2])
} else {
  stop("Provide start and end samples")
}
base_dir <- "/net/isi-scratch/kieran/"
load_all(paste0(base_dir, "switch/sctools/"))

# Load pseudotime assignments ------
if(!require(rhdf5)) stop("Need some hdf5 love")

source(paste0(base_dir, "GP/pseudogp2/stan/gbio/varygamma/diffexpr/common.R"))

##--------------- Edit here
tracefile <- paste0(base_dir, "GP/pseudogp2/data/varygamma_traces.h5")

h5outputfile <- paste0(base_dir, "GP/pseudogp2/data/varygamma/g2.hdf5")
##--------------- End edit


pst <- h5read(tracefile, "g2/pst")

##------- MONOCLE ONLY
source(paste0(base_dir, "GP/pseudogp2/stan/diffexpr/monocle/prep_data.R"))
sce <- load_data(base_dir)
## end


## sanity checking
stopifnot(start > 0 && end > 0 && start < nrow(pst) && end < nrow(pst) && start < end)

doDiffExprAnalysis(sce, pst, start, end, h5outputfile)






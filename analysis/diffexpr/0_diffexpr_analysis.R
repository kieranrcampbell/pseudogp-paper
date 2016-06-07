
# Differential expression comparison for uncertainty in pseudotime
# kieran.campbell@sjc.ox.ac.uk

library(devtools)
library(scater)
library(embeddr)
library(switchde)

source("analysis/diffexpr/common.R")

start <- NULL
end <- NULL
tracefile <- h5outputfile <- sce_file <- statusfile <- NULL

args <- commandArgs(trailingOnly = TRUE)

if(length(args) > 0) {
  tracefile <- args[1]
  sce_file <- args[2]
  h5outputfile <- args[3]
  startend <- strsplit(args[4], "_")[[1]][2]
  startend <- strsplit(startend, "x")[[1]]
  start <- as.numeric(startend[1])
  end <- as.numeric(startend[2])
  statusfile <- args[5]
} else {
  stop("Provide start and end samples")
}

# Load pseudotime assignments ------
if(!require(rhdf5)) stop("Need some hdf5 love")

print(paste("Performing differential expression analysis with start end", start, end))

pst <- h5read(tracefile, "pst")

load(sce_file)
sce@lowerDetectionLimit <- 0.1

## sanity checking
stopifnot(start > 0 && end > 0 && start < nrow(pst) && end < nrow(pst) && start < end)

doDiffExprAnalysis(sce, pst, start, end, h5outputfile)

conn <- file(statusfile)
writeLines("Finished", conn)
close(conn)


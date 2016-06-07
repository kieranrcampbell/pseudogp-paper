
# Differential expression comparison for uncertainty in pseudotime
# kieran.campbell@sjc.ox.ac.uk

library(devtools)
library(scater)
library(embeddr)
library(switchde)

source("analysis/diffexpr/common.R")

start <- NULL
end <- NULL
tracefile <- csv_file <- sce_file <- statusfile <- NULL

args <- commandArgs(trailingOnly = TRUE)

if(length(args) > 0) {
  tracefile <- args[1]
  sce_file <- args[2]
  csv_file <- args[3]
  startend <- strsplit(csv_file, "_")[[1]][2]
  startend <- strsplit(startend, ".csv")[[1]][1]
  startend <- strsplit(startend, "x")[[1]]
  start <- as.numeric(startend[1])
  end <- as.numeric(startend[2])
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

doDiffExprAnalysis(sce, pst, start, end, csv_file)

conn <- file(statusfile)
writeLines("Finished", conn)
close(conn)


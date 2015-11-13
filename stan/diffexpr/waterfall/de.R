
# Differential expression comparison for uncertainty in pseudotime
# kieran.campbell@sjc.ox.ac.uk

library(devtools)
library(scater)
library(embeddr)

doDiffExprAnalysis <- function(sce, pst, start, end, h5outputfile) {
  ## subset to genes expressed in at lesat 10% of cells
  n_cells_exprs <- rowSums(exprs(sce) > sce@lowerDetectionLimit)
  genes_to_use <- n_cells_exprs > (0.1 * ncol(sce)) # select genes expressed in at least 10% of cells
  sce <- sce[genes_to_use,]
  
  # What range of pst are we looking at?
  pst <- pst[start:end,] 
  
  ngenes <- nrow(sce) # number of genes
  n_pst <- nrow(pst) # number of pseudotime samples
  
  ss_pvals <- switch_pvals <- matrix(NA, nrow = ngenes, ncol = n_pst)

  for(i in 1:n_pst) {
    sce$pseudotime <- t_i <- pst[i,]
    
    ss_test <- pseudotime_test(sce, n_cores = 1)
    ss_pvals[,i] <- ss_test$p_val
    
    # switch model
    switch_pv <- testDE(sce)
    switch_pvals[,i] <- switch_pv[1,]
  }

  ## write results
  if(!file.exists(h5outputfile)) h5createFile(h5outputfile)
  h5createGroup(h5outputfile, "ss")
  h5createGroup(h5outputfile, "switch")
  h5write(ss_pvals, h5outputfile, paste0("ss/", as.character(start), "_", as.character(end)))
  h5write(switch_pvals, h5outputfile, paste0("switch/", as.character(start), "_", as.character(end)))
  print("[de.R] Done")
}

start <- NULL
end <- NULL

args <- commandArgs(trailingOnly = TRUE)

if(length(args) > 0) {
  start <- as.numeric(args[1])
  end <- as.numeric(args[2])
} else {
  stop("Provide start and end samples")
}
rdir <- "/net/isi-scratch/kieran/"
load_all(paste0(rdir, "switch/sctools/"))

# Load pseudotime assignments ------
if(!require(rhdf5)) stop("Need some hdf5 love")



##--------------- Edit here
tracefile <- paste0(rdir, "GP/pseudogp2/data/waterfall_stan_traces.h5")
scefile <- paste0(rdir, "datasets/waterfall/waterfall.Rdata")
h5outputfile <- paste0(rdir, "GP/pseudogp2/data/waterfall_diffexpr.h5")
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






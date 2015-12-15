
library(readr)
library(ggplot2)
library(cowplot)
library(MCMCglmm)
library(coda)
library(matrixStats)
library(reshape2)
library(scater)
library(embeddr)
library(dplyr)
library(rhdf5)

pvalsFromHDF5 <- function(h5file) {
  ssp <- h5read(h5file, "ss")
  swp <- h5read(h5file, "switch")
  
  ss_p <- do.call("cbind", ssp)
  switch_p <- do.call("cbind", swp)
  
  
  ## normal_switch gives a p-value of -1 if the alternative model couldn't be fit
  ## and -1 if the null model couldn't be fit
  optfail_per_gene <- apply(switch_p, 1, function(x) sum(x == -1))
  nullfail_per_gene <- apply(switch_p, 1, function(x) sum(x == -2)) ## this should be none
  
  stopifnot(sum(nullfail_per_gene) == 0) # we're in trouble if any of the null models failed
  
  switch_p_copy <- switch_p
  switch_p[switch_p == -1] <- 1
  
  ss_q <- apply(ss_p, 2, p.adjust, "BH")
  switch_q <- apply(switch_p, 2, p.adjust, "BH")
  
  ss_is_sig <- ss_q < 0.05
  switch_is_sig <- switch_q < 0.05
  
  ss_med_q <- rowMedians(ss_q)
  sw_med_q <- rowMedians(switch_q)
  
  ss_prop_sig <- rowSums(ss_is_sig) / ncol(ss_p)
  switch_prop_sig <- rowSums(switch_is_sig) / ncol(switch_p)
  
  return(list(sw_p = switch_p, ss_p = ss_p,
              sw_q = switch_q, ss_q = ss_q,
              sw_is_sig = switch_is_sig, ss_is_sig = ss_is_sig,
              sw_med_q = sw_med_q, ss_med_q = ss_med_q,
              sw_prop_sig = switch_prop_sig, ss_prop_sig = ss_prop_sig))
}


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

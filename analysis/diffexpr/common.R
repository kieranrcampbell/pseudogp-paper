
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
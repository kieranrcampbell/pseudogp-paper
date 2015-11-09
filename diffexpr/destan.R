
# Differential expression comparison for uncertainty in pseudotime
# kieran.campbell@sjc.ox.ac.uk

library(ggplot2)
library(readr)
library(devtools)
library(doParallel)
library(foreach)
library(scater)
library(embeddr)


rdir <- "/net/isi-scratch/kieran/"
load_all(paste0(rdir, "switch/sctools/"))
# load_all(paste0(rdir, "embeddr/embeddr"))
# load_all(paste0(rdir, "scater/scater"))

cluster <- to_keep <- pseudotimes <- NULL

# Load pseudotime assignments ------
if(require(rhdf5)) { # woo hdf5 is actually installed
  h5file <- paste0(rdir, "GP/gpseudotime/data/pseudotime_traces_6e6iterations_5e5gamma_5_10_2015.h5")
  tracefile <- paste0(rdir, "GP/pseudogp2/data/stan_traces.h5")
  to_keep <- h5read(h5file, "to_keep")
  pst <- pseudotimes <- h5read(tracefile, "pst")
  cluster <- h5read(h5file, "embeddr_cluster")
} else {
  stop("implement this properly")
  pseudotime_file <- paste0(rdir, "GP/gpseudotime/data/decsv/pseudotimes.csv")
  to_keep_file <- paste0(rdir, "GP/gpseudotime/data/decsv/to_keep.csv")
  cluster_file <- paste0(rdir, "GP/gpseudotime/data/decsv/cluster.csv")
  
  pseudotimes <- as.matrix(read_csv(pseudotime_file))
  to_keep <- read_csv(to_keep_file)$to_keep
  cluster  <- read_csv(cluster_file)$cluster
}

# Load gene expression data ---------
sce <- NULL
hsmm_data_available <- data(package='HSMMSingleCell')$results[,3]
if("HSMM" %in% hsmm_data_available) {
  data(HSMM)
  sce <- fromCellDataSet(HSMM, use_exprs_as = 'fpkm')
} else if("HSMM_expr_matrix" %in% hsmm_data_available) {
  library(HSMMSingleCell)
  data(HSMM_expr_matrix)
  data(HSMM_gene_annotation)
  data(HSMM_sample_sheet)
  
  pd <- new('AnnotatedDataFrame', data = HSMM_sample_sheet)
  fd <- new('AnnotatedDataFrame', data = HSMM_gene_annotation)
  sce <- newSCESet(fpkmData = HSMM_expr_matrix, phenoData = pd, featureData = fd)
} else {
  stop('No recognised data types in HSMMSingleCell')
}

sce@lowerDetectionLimit <- 0.1

# Subset off to the ones we want
sce <- sce[, cluster %in% 1:2]
sce <- sce[, which(to_keep == 1)]
sce <- calculateQCMetrics(sce)

# Look at what genes we want in our differential expression analysis -----

n_cells_exprs <- rowSums(exprs(sce) > sce@lowerDetectionLimit)

# Look at what genes we want in our differential expression analysis -----
# qplot(fData(sce)$n_cells_exprs)
genes_to_use <- n_cells_exprs > (0.1 * ncol(sce)) # select genes expressed in at least 10% of cells
sce <- sce[genes_to_use,]

# Time for some differential expression ------
#pst <- pseudotimes[500:nrow(pseudotimes),]


# What range of pst are we looking at?
args <- commandArgs(trailingOnly = TRUE)

start <- 1
end <- 2# nrow(pst)

if(length(args) > 0) {
  start <- as.numeric(args[1])
  end <- as.numeric(args[2])
  
  ## sanity checking
  stopifnot(start > 0 && end > 0 && start < nrow(pst) && end < nrow(pst) && start < end)
}

pst <- pst[start:end,] 

ngenes <- nrow(sce)
n_pst <- nrow(pst)

ss_pvals <- switch_pvals <- matrix(NA, nrow = ngenes, ncol = n_pst)

cl <- makeCluster(4)
registerDoParallel(cl)

ii <- 0


# print("do we get this far?")

# foreach(i = 1:n_pst) %dopar% {
for(i in 1:n_pst) {
  ii <- i
  sce$pseudotime <- t_i <- pst[i,]
  
  ss_test <- pseudotime_test(sce, n_cores = 1)
  ss_pvals[,i] <- ss_test$p_val
  
  # switch model
  switch_pv <- testDE(sce)
  switch_pvals[,i] <- switch_pv[1,]
}

if(require(rhdf5)) {
  h5outputfile <- paste0(rdir, "GP/pseudogp2/data/stan_diffexpr.h5")
  #h5createGroup(h5outputfile, "ss")
  #h5createGroup(h5outputfile, "switch")
  h5write(ss_pvals, h5outputfile, paste0("ss/", as.character(start), "_", as.character(end)))
  h5write(switch_pvals, h5outputfile, paste0("switch/", as.character(start), "_", as.character(end)))
} else {
  stop("reimplement!")
  ss_pval_output <- paste0(rdir, "GP/pseudogp2/data/ldl01/ss_pval_output", as.character(start), as.character(end), ".csv")
  switch_pval_output <- paste0(rdir, "GP/pseudogp2/data/ldl01/switch_pval_output", as.character(start), as.character(end), ".csv")
  
  write_csv(data.frame(ss_pvals), ss_pval_output)
  write_csv(data.frame(switch_pvals), switch_pval_output)
}








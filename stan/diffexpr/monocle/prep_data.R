
## Loads sce and pseudotimes and filters on genes etc

library(scater)
library(monocle)

load_data <- function(rdir = "/net/isi-scratch/kieran/") {
  cluster <- to_keep <- pseudotimes <- NULL
  
  # Load pseudotime assignments ------
  if(require(rhdf5)) { # woo hdf5 is actually installed
    h5file <- paste0(rdir, "GP/gpseudotime/data/pseudotime_traces_6e6iterations_5e5gamma_5_10_2015.h5")
    h5ls(h5file)
    
    to_keep <- h5read(h5file, "to_keep")
    pseudotimes <- h5read(h5file, "tchain")
    cluster <- h5read(h5file, "embeddr_cluster")
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
    sce <- newSCESet(fpkmData = HSMM_expr_matrix, phenoData = pd, 
                     featureData = fd, lowerDetectionLimit = 0.1)
  } else {
    stop('No recognised data types in HSMMSingleCell')
  }

  sce@lowerDetectionLimit <- 0.1
  
  # Subset off to the ones we want
  sce <- sce[, cluster %in% 1:2]
  sce <- sce[, which(to_keep == 1)]
  sce <- calculateQCMetrics(sce)
  
  n_cells_exprs <- rowSums(exprs(sce) > sce@lowerDetectionLimit)
  
  # Look at what genes we want in our differential expression analysis -----
  # qplot(fData(sce)$n_cells_exprs)
  genes_to_use <- n_cells_exprs > (0.1 * ncol(sce)) # select genes expressed in at least 10% of cells
  sce <- sce[genes_to_use,]

  return( sce )
}
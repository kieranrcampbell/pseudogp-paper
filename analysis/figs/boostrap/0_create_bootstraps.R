#' This script produces reduced dimension representations of 
#' the Trapnell et al dataset. Currently uses PCA with 80%
#' of cells subsampled.
#' 
#' We also save a full-cell representation too.

library(scater)
library(princurve)
library(dplyr)
library(readr)

set.seed(123L)

# We'll use the trapnell dataset for this
load("data/sce_trapnell.Rdata")
sce@featureControlInfo <- AnnotatedDataFrame()

ncells <- ncol(sce)

n_resamples <- 500

PCA_reps <- lapply(seq_len(n_resamples), function(i) {
  which_cells <- sample(ncells, ncells, replace = TRUE)
  new_exprs <- exprs(sce)[, which_cells]
  colnames(new_exprs) <- paste0("Cell", 1:ncol(new_exprs))
  new_pdata <- pData(sce)[which_cells,]
  rownames(new_pdata) <- colnames(new_exprs)
  new_sceset <- newSCESet(exprsData = new_exprs, phenoData = AnnotatedDataFrame(new_pdata))
  pca <- redDim(plotPCA(new_sceset, return_SCESet = TRUE))[,1:2]
  return(list(which_cells = which_cells, pca = pca))
})


# Now fit pseudotimes -----------------------------------------------------

principal_curve_fits <- sapply(PCA_reps, function(pcr) principal.curve(pcr$pca)$lambda)

which_cells <- sapply(PCA_reps, `[[`, "which_cells") %>% as_data_frame()
names(which_cells) <- paste0("Bootstrap", 1:ncol(which_cells))

pcf <- as_data_frame(principal_curve_fits)
names(pcf) <- names(which_cells)

write_csv(which_cells, "data/bootstrap/which_cells.csv")
write_csv(pcf, "data/bootstrap/bootstrapped_pseudotimes.csv")


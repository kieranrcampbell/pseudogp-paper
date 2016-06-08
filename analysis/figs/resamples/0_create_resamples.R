library(scater)

set.seed(123L)

# We'll use the trapnell dataset for this
load("data/sce_trapnell.Rdata")

ncells <- ncol(sce)
pct_subset <- 0.8
to_sample <- round(pct_subset * ncells)
n_resamples <- 100

PCA_reps <- lapply(seq_len(n_resamples), function(i) {
  which_cells <- sample(ncells, to_sample)
  pca <- redDim(plotPCA(sce[, which_cells], return_SCESet = TRUE))[,1:2]
  return(list(which_cells = which_cells, pca = pca))
})

save(PCA_reps, file = "data/resamples/pca_resamples.Rdata")

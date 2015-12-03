library(monocle)
library(scater)
library(rhdf5)

base_dir <- "~/mount/"
source(paste0(base_dir, "GP/pseudogp2/stan/diffexpr/monocle/prep_data.R"))
sce <- load_data(base_dir)

h5file <- paste0(base_dir, "GP/pseudogp2/data/5m_run_with_tau_traces.h5")
X <- h5read(h5file, "X")

sce <- plotPCA(sce, return_SCESet = TRUE)
plotReducedDim(sce)

cds <- toCellDataSet(sce)
sce2 <- fromCellDataSet(cds, use_exprs_as = "exprs")
plotReducedDim(sce2)

set.seed(123)
nsamples <- 20
standardize <- function(t) (t - min(t))/(max(t) - min(t))
ncells <- ncol(sce)
pseudosamples <- vector("list", ncells)
gt <- orderCells(cds)$Pseudotime
gt <- standardize(gt)

# did you ever see such awful R in your life
for(i in 1:nsamples) {
  chosen <- sample(ncells, round(0.6*ncells), replace=FALSE)
  cds_copy <- cds[,chosen]
  cds_copy@reducedDimS <- cds@reducedDimS[,chosen]
  cds_copy <- orderCells(cds_copy)
  pseudotimes <- standardize(cds_copy$Pseudotime)
  ## need to work out if we need to reverse the pseudotimes
  if(cor(gt[chosen], pseudotimes) < 0) pseudotimes <- 1 - pseudotimes
  
  for(j in 1:length(chosen)) {
    pseudosamples[[ chosen[j] ]] <- c(pseudosamples[[ chosen[j] ]], pseudotimes[j])
  }
}

pst_means <- sapply(pseudosamples, mean)
pst_sd <- sapply(pseudosamples, sd)
pst_95 <- 2*pst_sd

plot(gt, pst_means)
qplot(pst_95, geom='density')





#' We want to subset the SCESet to a computationally feasible number of cells
#' 

library(scater)

set.seed(123)

load("data/sce_trapnell.Rdata")

# First subset down to genes expressed in at least 10% of cells with
# minimum mean expression of 0.1:

prop_cells_exprs <- rowMeans(exprs(sce) > 0)

to_use <- prop_cells_exprs > 0.1 & rowMeans(exprs(sce)) > 0.1

sce <- sce[which(to_use), ]

ngenes <- nrow(sce)

to_de <- min(5000, ngenes)

to_use <- sample(ngenes, to_de)

sce <- sce[to_use, ]

save(sce, file = "data/resamples/sce_trapnell_resamples.Rdata")
#' This script produces reduced dimension representations of 
#' the Trapnell et al dataset. Currently uses PCA with 80%
#' of cells subsampled.
#' 
#' We also save a full-cell representation too.

library(scater)
library(dplyr)
library(readr)
library(pseudogp)

set.seed(123L)

# We'll use the trapnell dataset for this
load("data/sce_trapnell.Rdata")
sce@featureControlInfo <- AnnotatedDataFrame()

sce <- plotPCA(sce, return_SCESet = TRUE)

fit <- fitPseudotime(redDim(sce)[,1:2], "pca", iter = 5000, thin = 5,
                     smoothing_alpha = 12, smoothing_beta = 2)

pdf("figs/bootstrap_diagnostics.pdf")
plotDiagnostic(fit)
posteriorCurvePlot(redDim(sce)[,1:2], fit)
dev.off()

pseudotime_traces <- rstan::extract(fit, "t")$t %>% t() %>% tbl_df()

write_csv(pseudotime_traces, "data/bootstrap/gplvm_pseudotimes.csv")


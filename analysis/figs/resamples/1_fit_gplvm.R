#' This function fits GPLVM pseudotime curves to each
#' subsampled reduced dimension representation of the cells.

library(pseudogp)
library(methods)

set.seed(123L)

args <- commandArgs(trailingOnly = TRUE)

i <- as.numeric(args[1]) # which resample are we performing inference on?

load("data/resamples/pca_resamples.Rdata")

rs <- PCA_reps[[i]]
X <- rs$pca


fit <- fitPseudotime(X, "pca", iter = 5000, thin = 5,
smoothing_alpha = 10, smoothing_beta = 2)


pdf(paste0("data/resamples/diagnostic_plots/diagnostic_", i, ".pdf"), width = 7)
plotDiagnostic(fit)
posteriorCurvePlot(X, fit)
dev.off()

save(fit, file = paste0("data/resamples/gplvm_fits/fit_", i, ".Rdata"))
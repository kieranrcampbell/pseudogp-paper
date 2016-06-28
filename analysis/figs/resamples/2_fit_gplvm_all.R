#' This function fits GPLVM pseudotime curves to each
#' subsampled reduced dimension representation of the cells.

library(pseudogp)
library(methods)

set.seed(123L)

load("data/resamples/pca_all.Rdata")

X <- pca_all


fit <- fitPseudotime(X, "pca", iter = 5000, thin = 5,
                    smoothing_alpha = 16, smoothing_beta = 2)


pdf("data/resamples/diagnostic_plots/diagnostic_all.pdf", width = 7)
plotDiagnostic(fit)
posteriorCurvePlot(X, fit)
dev.off()

save(fit, file = "data/resamples/gplvm_fit_all.Rdata")

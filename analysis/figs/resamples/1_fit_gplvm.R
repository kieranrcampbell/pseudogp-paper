library(pseudogp)

args <- commandArgs(trailingOnly = TRUE)

i <- as.numeric(args[1]) # which resample are we performing inference on?

load("data/resamples/pca_resamples.Rdata")

rs <- PCA_reps[[i]]
X <- rs$pca

fit <- fitPseudotime(X, "principal_curve", iter = 5000, thin = 5)

pdf(paste0("data/resamples/diagnostic_plots/diagnostic_", i, ".pdf"), width = 7)
plotDiagnostic(fit)
posteriorCurvePlot(X, fit)
dev.off()

save(fit, file = paste0("data/resamples/gplvm_fits/fit_", i, ".Rdata"))
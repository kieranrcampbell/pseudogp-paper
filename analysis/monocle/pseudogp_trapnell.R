#' ## Fitting probabilistic trajectories to Burns et al dataset
#' ### Kieran Campbell <kieran.campbell@sjc.ox.ac.uk>
#' 
#' This document goes through fitting the probabilistic curve to the Burns et al.
#' dataset. To generate this document in markdown, run `knitr::spin("pseudogp_trapnell.R")` 
#' from
#' within R.

#+ setup
library(rstan)
library(rhdf5)
library(MCMCglmm)
library(coda)
library(ggplot2)

devtools::load_all("~/oxford/pseudogp") # TODO - install and replace with `library(pseudogp)`
set.seed(123)


#' Now read in the data
#' 
#+ data-read
base_dir <- "~/mount"

h5file <- file.path(base_dir, "GP/pseudogp2/data/trapnell_embeddings.h5")
output_hdf5 <- file.path(base_dir, "GP/pseudogp2/data/monocle_stan_traces.h5")

X <- h5read(h5file, "Xle")
t_gt <- h5read(h5file, "t_gt")

#' Now fit the probabilistic pseudotime:
#+ fit-pseudotime, message = FALSE
fit <- fitPseudotime(X, smoothing_alpha = 30, smoothing_beta = 5, seed = 123)

#' Plot posterior mean curves
#+ posteriorcplt, fig.width=6, fig.height=5
posteriorCurvePlot(X, fit, posterior_mean = TRUE)

#' Posterior boxplot of pseudotime distribution
#+ posteriorbplt, fig.width=7, fig.height=5
posteriorBoxplot(fit)

#' Because we're dealing with stan objects we can very easily extract the posterior samples:
#+ extract samples
pst <- extract(fit, "t")
tmcmc <- mcmc(pst$t)

smcmc <- mcmc(extract(fit, "sigma")$sigma[,1,]) # need to slice middle index to get single rep.
lmcmc <- mcmc(extract(fit, "lambda")$lambda[,1,])

#' ...and save them to HDF5
#+ save-hdf5
if(!file.exists(output_hdf5)) h5createFile(output_hdf5)

h5write(X, output_hdf5, "X")
h5write(pst$t, output_hdf5, "pst")
h5write(as.matrix(smcmc), output_hdf5, "sigma")
h5write(as.matrix(lmcmc), output_hdf5, "lambda")
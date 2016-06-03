#' # Probabilsitic trajectory fitting for the Shin et al. (2015) dataset
#' ### Kieran Campbell <kieran.campbell@sjc.ox.ac.uk>
#' 
#' To turn this into markdown, run `knitr::spin("pseudogp_shin.R")` from
#' within R.
#' 
#+setup
library(rstan)
library(rhdf5)
library(MCMCglmm)
library(coda)
library(ggplot2)
library(moments)

base_dir <- "~/mount/"

h5file <- file.path(base_dir, "pseudogp-paper/data/waterfall_embeddings.h5")
output_hdf5 <- file.path(base_dir, "pseudogp-paper/data/waterfall_stan_traces.h5")

devtools::load_all("~/oxford/pseudogp") # TODO: change to library(pseudogp)
set.seed(123)


#' Read in the data
#+ read-data
X <- h5read(h5file, "Xpca")
wpst <- h5read(h5file, "wpst") # the pseudotime assigned using Waterfall

#' Fit the pseudotime model
#+ fit-pseudotime
fit <- fitPseudotime(X, initialise_from = "pca", smoothing_alpha = 8, smoothing_beta = 2, seed = 123)

#' Plot the posterior mean curve
#+ posmean-curve, fig.width=6, fig.height=5
posteriorCurvePlot(X, fit)

#' And the boxplots
#+ posmean-boxplot, fig.width=6, fig.height=5
posteriorBoxplot(fit)

#' We can extract the posterior traces and compare the MAP estimate to the waterfall fit
#+ compare-map, fig.width=5, fig.height=4
pst <- extract(fit, "t")
tmcmc <- mcmc(pst$t)
post_mean <- posterior.mode(tmcmc)
qplot(wpst, post_mean) + theme_bw() + xlab("Waterfall fit") + ylab("Pseudogp map")

smcmc <- mcmc(extract(fit, "sigma")$sigma[,1,])
lmcmc <- mcmc(extract(fit, "lambda")$lambda[,1,])

#' And finally save everything to HDF5
h5createFile(output_hdf5)

h5write(X, output_hdf5, "X")
h5write(pst$t, output_hdf5, "pst")
h5write(as.matrix(smcmc), output_hdf5, "sigma")
h5write(as.matrix(lmcmc), output_hdf5, "lambda")


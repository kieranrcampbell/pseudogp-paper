#' ## Fitting probabilistic trajectories to Burns et al dataset
#' ### Kieran Campbell <kieran.campbell@sjc.ox.ac.uk>
#' 
#' This document goes through fitting the probabilistic curve to the Burns et al.
#' dataset. To generate this document in markdown, run `knitr::spin("burns_pseudogp.R")` 
#' from
#' within R.


#' Quick setup
#' 
#+ setup
library(rstan)
library(rhdf5)
library(MCMCglmm)
library(coda)
library(ggplot2)
library(moments)
library(devtools)

devtools::load_all("~/oxford/pseudogp") # TODO - install and replace with `library(pseudogp)`
set.seed(123)

#' Now read in the data
#' 
#+ data-read
base_dir <- "~/mount"

h5file <- file.path(base_dir, "GP/pseudogp2/data/ear_embeddings.h5")
output_hdf5 <- file.path(base_dir, "GP/pseudogp2/data/ear_stan_traces.h5")

X <- h5read(h5file, "Xle")
X <- apply(X, 2, function(x) (x - mean(x)) / sd(x))
t_gt <- h5read(h5file, "t_gt")

#' and quickly plot to make sure it looks right
#+ quick-plot, fig.width=5, fig.height = 5
ggplot(data.frame(X, t_gt)) + 
  geom_point(aes(x = X1, y = X2, color = t_gt)) + theme_bw()

#' ### Fit the pseudotime
#+ pseudo-fit
fit <- fitPseudotime(X, initialise_from = "principal_curve", 
                     smoothing_alpha = 9, smoothing_beta = 1, seed = 123)

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




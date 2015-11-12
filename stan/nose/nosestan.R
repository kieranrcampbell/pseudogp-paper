library(rstan)
library(rhdf5)
library(MCMCglmm)
library(coda)
library(ggplot2)
library(moments)

set.seed(123)

setwd("/net/isi-scratch/kieran/GP/pseudogp2/")
source("gputils//gputils.R")
h5file = "data/nose_embeddings.h5"
X <- h5read(h5file, "X")
X <- apply(X, 2, function(x) (x - mean(x)) / sd(x))
wpst <- h5read(h5file, "wpst")


data <- list(X = X, N = nrow(X))

fit <- stan(file = "stan/nose/pseudogp_nose.stan", data = data, 
            iter = 1000, chains = 1)

plot(fit, pars = "t")
plot(fit, pars = "lambda")
plot(fit, pars = "g")

pst <- extract(fit, "t")

tmcmc <- mcmc(pst$t)
post_mean <- posterior.mode(tmcmc)
plot(wpst, post_mean)

smcmc <- mcmc(extract(fit, "sigma")$sigma)
lmcmc <- mcmc(extract(fit, "lambda")$lambda)

smap <- posterior.mode(smcmc)
lmap <- posterior.mode(lmcmc)

plot_posterior_mean(X, wpst, post_mean, lmap, smap)

h5file <- "/net/isi-scratch/kieran/GP/pseudogp2/data/nose_stan_traces.h5"
if(!file.exists(h5file)) h5createFile(h5file)

h5write(X, h5file, "X")
h5write(pst$t, h5file, "pst")
h5write(as.matrix(smcmc), h5file, "sigma")
h5write(as.matrix(lmcmc), h5file, "lambda")


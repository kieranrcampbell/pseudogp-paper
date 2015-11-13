library(rstan)
library(rhdf5)
library(MCMCglmm)
library(coda)
library(ggplot2)
library(moments)

set.seed(123)

setwd("~/mount/GP/pseudogp2/")
source("gputils//gputils.R")
h5file = "data/waterfall_embeddings.h5"
X <- h5read(h5file, "X")
X <- apply(X, 2, function(x) (x - mean(x)) / sd(x))
wpst <- h5read(h5file, "wpst")


data <- list(X = X, N = nrow(X))

ti <- X[,1]
ti <- (ti - min(ti) + 10e-4) / (max(ti) - min(ti) + 10e-3)
init <- list(t = ti, lambda = c(5,5), sigma = c(0.2, 0.2), g = 1)
l_init <- lapply(1:4, function(x) init)
fit <- stan(file = "stan/waterfall/pseudogp_waterfall.stan", data = data, 
            iter = 1000, chains = 4, init = l_init)

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

h5file <- "~/mount/GP/pseudogp2/data/waterfall_stan_traces.h5"
if(!file.exists(h5file)) h5createFile(h5file)

inds <- sample(nrow(tmcmc), 500)

h5write(X, h5file, "X")
h5write(pst$t[inds,], h5file, "pst")
h5write(as.matrix(smcmc)[inds,], h5file, "sigma")
h5write(as.matrix(lmcmc)[inds,], h5file, "lambda")


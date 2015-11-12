library(rstan)
library(rhdf5)
library(MCMCglmm)
library(coda)
library(ggplot2)
library(moments)

set.seed(123)

setwd("/net/isi-scratch/kieran/GP/pseudogp2/")
source("gputils//gputils.R")
h5file = "data/ear_embeddings.h5"
X <- h5read(h5file, "X")
X <- apply(X, 2, function(x) (x - mean(x)) / sd(x))
t_gt <- h5read(h5file, "t_gt")


Xp <- X
Xp <- X[X[,2] < 2,]
Xp <- apply(Xp, 2, function(x) (x - mean(x)) / sd(x))

data <- list(X = Xp, N = nrow(Xp))

fit <- stan(file = "stan/ear//pseudogp_ear.stan", data = data, 
            iter = 1000, chains = 1)

plot(fit, pars = "t")
plot(fit, pars = "lambda")
plot(fit, pars = "g")

pst <- extract(fit, "t")

tmcmc <- mcmc(pst$t)
post_mean <- posterior.mode(tmcmc)
plot(t_gt[X[,2] < 2], post_mean)

smcmc <- mcmc(extract(fit, "sigma")$sigma)
lmcmc <- mcmc(extract(fit, "lambda")$lambda)

smap <- posterior.mode(smcmc)
lmap <- posterior.mode(lmcmc)

plot_posterior_mean(Xp, post_mean, post_mean, lmap, smap)

h5file <- "/net/isi-scratch/kieran/GP/pseudogp2/data/ear_stan_traces.h5"
if(!file.exists(h5file)) h5createFile(h5file)

h5write(Xp, h5file, "X")
h5write(pst$t, h5file, "pst")
h5write(as.matrix(smcmc), h5file, "sigma")
h5write(as.matrix(lmcmc), h5file, "lambda")
h5write(X[,2] < 2, h5file, "to_keep")




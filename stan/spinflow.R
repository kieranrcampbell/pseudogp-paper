
#' ## Supplementary analysis for pseudogp
#' 
#' In this document we compare some odds and ends for the pseudogp paper, including:
#' * Effect of unconstrained prior
#' * Effect of prior variance in unconstrained case
#' * Maximum likelihood estimate

#' First - setup:

#+ setup, cache=TRUE, message = FALSE
library(rstan)
library(rhdf5)
library(MCMCglmm)
library(coda)
library(ggplot2)
library(mvtnorm)
library(pscl)

set.seed(123)
setwd("/net/isi-scratch/kieran/GP/pseudogp2/stan")
source("../gputils//gputils.R")

h5file <- "/net/isi-scratch/kieran/GP/pseudogp2/data/5m_run_with_tau_traces.h5"

X <- h5read(h5file, "X")
X <- apply(X, 2, function(x) (x - mean(x)) / sd(x))
t_gt <- h5read(h5file, "t_gt")

#' We can plot the embedding:
#+ plot-embedding
ggplot(data.frame(X, t_gt)) + geom_point(aes(x = X1, y = X2, color = t_gt), size = 3) +
  scale_color_continuous(low = "darkred", high = "yellow") + theme_bw()

#' ### Constrained model

#+ constr-model, cache=TRUE, results='hide', warning = FALSE
data <- list(X = X, N = nrow(X))
fit <- stan(file = "pseudogp.stan", data = data, 
             iter = 1000, chains = 1)
plot(fit, pars = "t")

#' Plot model against principal curve fit
#+ pc-fit, message=FALSE, cache=TRUE
pst <- extract(fit, "t")$t
tmcmc <- mcmc(pst)
post_mean <- posterior.mode(tmcmc)
qplot(t_gt, post_mean, size = 3, alpha = 0.8) + xlab("Principal curve fit") + 
  ylab("Posterior mode") + theme_bw()


#' ### Maximum A-Posterior Estimates

#+ MLE, cache=TRUE
opt <- optimizing(get_stanmodel(fit), data = data)
t_opt <- opt$par[grep("t", names(opt$par))]
qplot(t_gt, t_opt, size = 3, alpha = 0.8) + xlab("Principal curve fit") +
  ylab("MAP estimate") + theme_bw()


#' Now modify the model slightly so the constrains on t are between 0 and 0.1:
#+ shrink-pst, message=FALSE, cache=TRUE, results='hide'
gamma_alpha <- 1 
gamma_beta <- 0.7
data <- list(X = X, N = nrow(X),
             t_lower = 0, t_upper = 0.1, 
             tmean = 0.05, tvar = 1,
             gamma_alpha = gamma_alpha, gamma_beta = gamma_beta)
shrink_fit <- stan(file = "models/pgp_adjustable.stan", data = data, iter = 1000, chains = 1)
plot(shrink_fit, pars = "t")
plot(shrink_fit, pars = "lambda")
plot(shrink_fit, pars = "g")
pst <- extract(shrink_fit, "t")$t
tmcmc <- mcmc(pst)
post_mean <- posterior.mode(tmcmc)
qplot(t_gt, post_mean, size = 3, alpha = 0.8) + xlab("Principal curve fit") + 
  ylab("Posterior mode") + theme_bw()

#' So picks up posterior fine - but you seriously need to decrease the strength of the
#' prior on lambda to get it to not just sample from the posterior. Look at MAP estimate:
#+ map-shrink, message=FALSE, cache=TRUE
opt <- optimizing(get_stanmodel(shrink_fit), data = data)
t_opt <- opt$par[grep("t", names(opt$par))]
qplot(t_gt, t_opt, size = 3, alpha = 0.8) + xlab("Principal curve fit") +
  ylab("MAP estimate") + theme_bw()

#' Evaluate the log-likelihood of the traces:
#+ evalloglik, cache=TRUE
eval_loglik <- function(X, t, lambda, sigma, gamma, data) {
  ll <- 0
  Sigma <- lapply(1:2, function(i) cov_matrix(t, t, lambda = lambda[i], sigma = sigma[i]))
  ll <- ll + sum(sapply(1:2, function(i) dmvnorm(X[,i], sigma = Sigma[[i]], log = TRUE)))
  ll <- ll + sum(log(sapply(sigma, densigamma, 1.0, 1.0)))
  ll <- ll + sum(sapply(lambda, dexp, gamma, TRUE))
  ll <- ll + dgamma(gamma, shape = data$gamma_alpha, scale = data$gamma_beta, log = TRUE)
  ll <- ll + sum(sapply(t, dnorm, data$tmean, data$tvar, TRUE))
  return(as.numeric(ll))
}

lambda <- opt$par[grep("lambda", names(opt$par))]
sigma <- opt$par[grep("sigma", names(opt$par))]
gamma <- opt$par[grep("^g$", names(opt$par))]
eval_loglik(X, t_opt, lambda, sigma, gamma, data)

#' ### Randomized data
#' We can randomise the data and see if we still pick up the ordering. Do
#' it twice to see what happens:

#+ randomized, cache=TRUE, warning = FALSE, message = FALSE, results='hide'
set.seed(123)
randomized_fits <- lapply(1:2, function(i) {
  Xr <- cbind(sample(X[,1]), sample(X[,2]))
  data <- list(X = Xr, N = nrow(X))
  fit <- stan(file = "pseudogp.stan", data = data, 
              iter = 1000, chains = 1, seed = 123)
  pst <- extract(fit, "t")$t
  tmcmc <- mcmc(pst)
  lmap <- posterior.mode(mcmc(extract(fit, "lambda")$lambda))
  smap <- posterior.mode(mcmc(extract(fit, "sigma")$sigma))
  post_mean <- posterior.mode(tmcmc)
  return(list(Xr = Xr, post_mean = post_mean, lmap = lmap, smap = smap, fit = fit))
})

#' And plot the results:
#+ plot-randomized, cache=TRUE, warning = FALSE, message = FALSE
source("../gputils//gputils.R")
plotlist <- lapply(randomized_fits, function(f) plot_posterior_mean(f$Xr, t_gt, f$post_mean, f$lmap, f$smap, curve_color = "red"))

cowplot::plot_grid(plotlist = plotlist, ncol = 1)

plot(randomized_fits[[1]]$fit, pars = "t")
plot(randomized_fits[[2]]$fit, pars = "t")

                   
#' ### Unconstrained prior with multiple chains
  
#+ unconst-mult, cache=TRUE, message = FALSE, warning = FALSE, results='hide'
data <- list(X = X, N = nrow(X), prior_var = 1)

ufit <- stan(file = "pseudogp_unconstr.stan", data = data, 
            iter = 1000, chains = 4)
plot(ufit, pars = "t")

#' Now see how the pseudotimes compare to the principal curve fit:


#+ comparestrands, cache=TRUE
pst <- extract(ufit, "t", permute = FALSE)
pmcmc_mode <- data.frame(apply(pst, 2, function(x) posterior.mode(mcmc(x))))
pmcmc_mode$pc <- t_gt
m <- reshape2::melt(pmcmc_mode, id.vars = "pc")
ggplot(m) + geom_point(aes(x = pc, y = value)) + facet_wrap(~ variable) +
  theme_bw()
ggplot(m) + geom_point(aes(x = pc, y = value, color = variable)) + theme_bw()


#' Now increase the prior variance to 10:
#+ unconst-mult-incvar, cache=TRUE, message = FALSE, warning = FALSE, results='hide'
data <- list(X = X, N = nrow(X), prior_var = 10)

uvfit <- stan(file = "pseudogp_unconstr.stan", data = data, 
             iter = 1000, chains = 4)
plot(uvfit, pars = "t")

#' Now see how the pseudotimes compare to the principal curve fit:


#+ comparestrands-incvar, cache=TRUE
pst <- extract(uvfit, "t", permute = FALSE)
pmcmc_mode <- data.frame(apply(pst, 2, function(x) posterior.mode(mcmc(x))))
pmcmc_mode$pc <- t_gt
m <- reshape2::melt(pmcmc_mode, id.vars = "pc")
ggplot(m) + geom_point(aes(x = pc, y = value)) + facet_wrap(~ variable) +
  theme_bw()
ggplot(m) + geom_point(aes(x = pc, y = value, color = variable)) + theme_bw()



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

set.seed(123)

setwd("/net/isi-scratch/kieran/GP/pseudogp2/stan")

h5file <- "/net/isi-scratch/kieran/GP/pseudogp2/data/5m_run_with_tau_traces.h5"

X <- h5read(h5file, "X")
X <- apply(X, 2, function(x) (x - mean(x)) / sd(x))
t_gt <- h5read(h5file, "t_gt")

#' We can plot the embedding:
#+ plot-embedding
ggplot(data.frame(X, t_gt)) + geom_point(aes(x = X1, y = X2, color = t_gt), size = 3) +
  scale_color_continuous(low = "darkred", high = "yellow") + theme_bw()

#' ### Constrained model

#+ constr-model, cache=TRUE
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

#' ### Unconstrained prior with multiple chains
  
#+ unconst-mult, cache=TRUE, message = FALSE, warning = FALSE
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
#+ unconst-mult-incvar, cache=TRUE, message = FALSE, warning = FALSE
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


#' ### Maximum A-Posterior Estimates

#+ MLE, cache=TRUE
opt <- optimizing(get_stanmodel(fit), data = data)
t_opt <- opt$par[grep("t", names(opt$par))]
qplot(t_gt, t_opt, size = 3, alpha = 0.8) + xlab("Principal curve fit") +
  ylab("MAP estimate") + theme_bw()

library(rstan)
library(rhdf5)
library(MCMCglmm)
library(coda)
library(ggplot2)
library(moments)

set.seed(123)

setwd("/net/isi-scratch/kieran/GP/pseudogp2/")
h5file = "data/ear_embeddings.h5"
X <- h5read(h5file, "X")
X <- apply(X, 2, function(x) (x - mean(x)) / sd(x))
t_gt <- h5read(h5file, "t_gt")

data <- list(X = X, N = nrow(X))

fit <- stan(file = "stan/pseudogp.stan", data = data, 
            iter = 2000, chains = 1)

plot(fit, pars = "t")
plot(fit, pars = "lambda")
plot(fit, pars = "g")

pst <- extract(fit, "t")

#tmcmc <- mcmc(pst$t)
tmcmc <- mcmc(pst$t)
post_mean <- posterior.mode(tmcmc)
plot(t_gt, post_mean)

smcmc <- mcmc(extract(fit, "sigma")$sigma)
lmcmc <- mcmc(extract(fit, "lambda")$lambda)

smap <- posterior.mode(smcmc)
lmap <- posterior.mode(lmcmc)

plot_posterior_mean(post_mean, lmap, smap)

cov_matrix <- function(t1, t2, lambda, sigma = NULL) {
  n1 <- length(t1)
  n2 <- length(t2)
  C <- matrix(NA, nrow = n1, ncol = n2)
  for(i in 1:n1) {
    for(j in 1:n2) {
      C[i, j] <- exp(-lambda * (t1[i] - t2[j])^2)
    }
  }  
  if(!is.null(sigma)) {
    stopifnot(n1 == n2)
    C <- C + sigma * diag(n1)
  }
  return ( C )
}

plot_posterior_mean <- function(t, l, s, nnt = 80) {
  nt <- runif(nnt)
  K_y <- lapply(1:2, function(i) cov_matrix(t, t, as.numeric(l[i]), as.numeric(s[i])))
  K_star <- lapply(1:2, function(i) cov_matrix(t, nt, as.numeric(l[i])))
  K_dstar <- lapply(1:2, function(i) cov_matrix(nt, nt, as.numeric(l[i])))
  
  mu_star <- lapply(1:2, function(i) {
    t(K_star[[i]]) %*% solve(K_y[[i]]) %*% X[,i]
  })
  
  mus <- do.call(cbind, mu_star)
  pdf <- data.frame(mus[order(nt),], nt = nt[order(nt)])
  ggplot() + 
    geom_point(data = data.frame(X, t_gt), aes(x = X1, y = X2, color = t_gt), size = 3, alpha = 0.5) + 
    geom_path(data = pdf, aes(x = X1, y = X2, color = nt), size = 2, alpha = .8) + theme_bw() +
    geom_rug(data = data.frame(X, t_gt), aes(x = X1, y = X2))
}


# Reconstruct covariance matrix -------------------------------------------

t <- post_mean ; l <- lmap; s <- smap
K_y <- lapply(1:2, function(i) cov_matrix(t, t, as.numeric(l[i]), as.numeric(s[i])))




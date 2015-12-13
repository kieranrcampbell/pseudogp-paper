library(rstan)
library(rhdf5)
library(MCMCglmm)
library(coda)
library(ggplot2)


base_dir <- "/net/isi-scratch/kieran/"
setwd(file.path(base_dir, "GP/pseudogp2/stan"))

h5file <- file.path(base_dir, "GP/pseudogp2/data/5m_run_with_tau_traces.h5")
#h5file = "~/mount/GP/pseudogp2/data/ear_embeddings.h5"
X <- h5read(h5file, "X")
X <- apply(X, 2, function(x) (x - mean(x)) / sd(x))
t_gt <- h5read(h5file, "t_gt")

data <- list(X = X, N = nrow(X))

xpca <- prcomp(X)$x[,1]
t0 <- (xpca - min(xpca)) / (max(xpca) - min(xpca))

init <- list(list(t = t0),
             lambda = c(1,1),
             sigma = c(1,1))

fit <- stan(file = "monocle/pseudogp.stan", data = data,
            iter = 1000, chains = 1)

# opt <- optimizing(get_stanmodel(fit), data = data)
# t_opt <- opt$par[grep("t", names(opt$par))]
# plot(t_gt, t_opt)

plot(fit, pars = "t")
plot(fit, pars = "lambda")
plot(fit, pars = "g")

pst <- extract(fit, "t")

#tmcmc <- mcmc(pst$t)
tmcmc <- mcmc(pst$t)
post_mean <- posterior.mode(tmcmc)
plot(t_gt, post_mean)

rstan::traceplot(fit, pars = "lambda")
plot(tmcmc, pars = "t")

lmcmc <- mcmc(extract(fit, "lambda")$lambda)
smcmc <- mcmc(extract(fit, "sigma")$sigma)


## plot the predictive mean


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

t <- posterior.mode(tmcmc); l <- posterior.mode(lmcmc) ; s <- posterior.mode(smcmc)
plot_posterior_mean(t, l, s, nnt = 200)

ind <- sample(nrow(tmcmc), 4)
plts <- lapply(ind, function(i) plot_posterior_mean(tmcmc[i,], lmcmc[i,], smcmc[i,]))
cowplot::plot_grid(plotlist = plts)


tracefile <- "~/mount/GP/pseudogp2/data/stan_traces.h5"
h5createFile(tracefile)
h5write(pst$t, tracefile, "pst")

K_y <- lapply(1:2, function(i) cov_matrix(t, t, as.numeric(l[i]), as.numeric(s[i])))

library(rstan)
library(rhdf5)
library(MCMCglmm)
library(coda)
library(ggplot2)


base_dir <- "/net/isi-scratch/kieran/"
setwd(paste0(base_dir, "GP/pseudogp2/stan/multi"))
source("../../diffexpr/prep_data.R")


h5file <- paste0(base_dir, "GP/pseudogp2/data/5m_run_with_tau_traces.h5")
#h5file = "~/mount/GP/pseudogp2/data/ear_embeddings.h5"


# Laplacian eigenmaps representation --------------------------------------
X <- h5read(h5file, "X") 
X <- apply(X, 2, function(x) (x - mean(x)) / sd(x)) # this is the LE representation
t_gt <- h5read(h5file, "t_gt") # principal curves rep


# PCA representation ------------------------------------------------------
sce <- load_data()$sce
sce$pseudotime <- t_gt
sce <- plotPCA(sce, colour_by = "pseudotime", return_SCESet = TRUE)
Y <- redDim(sce)
Y <- apply(Y, 2, function(x) (x - mean(x)) / sd(x)) # our PCA embedding


# T-SNE representation ----------------------------------------------------

set.seed(1234)
sce <- plotTSNE(sce, colour_by = "pseudotime", perplexity = 3, return_SCESet = TRUE)
Z <- redDim(sce)
Z <- apply(Z, 2, function(x) (x - mean(x)) / sd(x)) # our t-SNE embedding
stopifnot(all(dim(X) == dim(Y)))
stopifnot(all(dim(Z) == dim(X)))

Ns <- 3
dx <- array(dim = c(Ns, 2, nrow(X)))
dx[1,,] <- t(X)
dx[2,,] <- t(Y)
dx[3,,] <- t(Z)

data <- list(Ns = Ns, P = 2, N = nrow(X), X = dx)


fit <- stan(file = "pgparray.stan", data = data, 
            iter = 1000, chains = 1)

plot(fit, pars = "t")
plot(fit, pars = "lambda")
plot(fit, pars = "g")

pst <- extract(fit, "t")

#tmcmc <- mcmc(pst$t)
tmcmc <- mcmc(pst$t)
post_mean <- posterior.mode(tmcmc)
plot(t_gt, post_mean)

lmcmc <- lapply(1:Ns, function(i) mcmc(extract(fit, "lambda")$lambda[,i,]))
lmap <- lapply(lmcmc, posterior.mode)

smcmc <- lapply(1:Ns, function(i) mcmc(extract(fit, "sigma")$sigma[,i,]))
smap <- lapply(smcmc, posterior.mode)


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

plot_posterior_mean <- function(X, t_gt, t, l, s, nnt = 80, reverse = FALSE) {
  if(reverse) t_gt <- 1 - t_gt
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
    geom_path(data = pdf, aes(x = X1, y = X2, color = nt), size = 2, alpha = .8) + theme_bw()
}

t <- posterior.mode(tmcmc)
colnames(X) <- colnames(Y) <- colnames(X)  <- NULL
xx <- list(X, Y, Z)

plots <- lapply(1:Ns, function(i) plot_posterior_mean(xx[[i]], t_gt, t, lmap[[i]], smap[[i]]))

joint_plot <- plot_grid(plotlist = plots, ncol = 1, labels = c("Laplacian eigenmaps", "PCA", "t-SNE"))
cowplot::ggsave("joint.png", joint_plot, width = 4, height = 3, scale = 2)

ind <- sample(nrow(tmcmc), 4)
plts <- lapply(ind, function(i) plot_posterior_mean(tmcmc[i,], lmcmc[i,], smcmc[i,]))
cowplot::plot_grid(plotlist = plts)


tracefile <- "~/mount/GP/pseudogp2/data/stan_traces.h5"
h5createFile(tracefile)
h5write(pst$t, tracefile, "pst")

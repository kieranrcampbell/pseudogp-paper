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

stopifnot(all(dim(X) == dim(Y)))
data <- list(X = X, Y = Y, N = nrow(X))

fit <- stan(file = "pgpmulti.stan", data = data, 
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

lXmcmc <- mcmc(extract(fit, "lambdaX")$lambdaX)
lYmcmc <- mcmc(extract(fit, "lambdaY")$lambdaY)

sXmcmc <- mcmc(extract(fit, "sigmaX")$sigmaX)
sYmcmc <- mcmc(extract(fit, "sigmaY")$sigmaY)


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

plot_posterior_mean <- function(X, t_gt, t, l, s, nnt = 80) {
  t_gt <- 1 - t_gt
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
lX <- posterior.mode(lXmcmc) ; sX <- posterior.mode(sXmcmc)
lY <- posterior.mode(lXmcmc) ; sY <- posterior.mode(sXmcmc)

colnames(Y) <- NULL
le_plt <- plot_posterior_mean(X, t_gt, t, lX, sX, nnt = 200)
pca_plt <- plot_posterior_mean(Y, t_gt, t, lY, sY, nnt = 200)

plot_grid(le_plt, pca_plt, ncol = 1, labels = c("Laplacian eigenmaps", "PCA"))

ind <- sample(nrow(tmcmc), 4)
plts <- lapply(ind, function(i) plot_posterior_mean(tmcmc[i,], lmcmc[i,], smcmc[i,]))
cowplot::plot_grid(plotlist = plts)


tracefile <- "~/mount/GP/pseudogp2/data/stan_traces.h5"
h5createFile(tracefile)
h5write(pst$t, tracefile, "pst")

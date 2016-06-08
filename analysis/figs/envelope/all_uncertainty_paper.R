library(coda)
library(MCMCglmm)
library(mvtnorm)
library(ggplot2)
library(rhdf5)
library(cowplot)
library(pseudogp)


makeUncertaintyPlot <- function(X, tmap, lmap, smap, xl, yl, nnt = 100, nx = 100, h = 2) {
  nt <- runif(nnt)
  
  K_y <- lapply(1:2, function(i) cov_matrix(tmap, tmap, as.numeric(lmap[i]), as.numeric(smap[i])))
  K_star <- lapply(1:2, function(i) cov_matrix(tmap, nt, as.numeric(lmap[i])))
  K_dstar <- lapply(1:2, function(i) cov_matrix(nt, nt, as.numeric(lmap[i])))
  
  mu_star <- lapply(1:2, function(i) {
    t(K_star[[i]]) %*% solve(K_y[[i]]) %*% X[,i]
  })
  
  Sigma <- lapply(1:2, function(i) K_dstar[[i]] - t(K_star[[i]]) %*% solve(K_y[[i]]) %*% K_star[[i]] + diag(nnt) * smap[i])
  
  Xgen <- sapply(1:2, function(i) rmvnorm(nx, mu_star[[i]], Sigma[[i]]))
  
  dxgen <- data.frame(Xgen)
  dx <- data.frame(X)
  dm <- data.frame(do.call("cbind", mu_star))
  dm$t <- nt
  dm <- dplyr::arrange(dm, t)
  
  ggplot() +
    stat_density2d(data = dxgen, aes(x = X1, y = X2, fill = ..density.., alpha = ..density..), 
                   geom = "tile", n = 200, contour=FALSE, h = h) +
    geom_point(data = dx, aes(x = X1, y = X2), shape = 21,
               fill = 'black', colour = 'white', size = 3, alpha = 0.6) +
    theme(legend.position="none") + scale_fill_gradient (low = "#FFFFFF", high = "darkred") +
    geom_path(data = dm, aes(x = X1, y = X2), size = 1.4, alpha = 0.7) +
    xlab("Component 1") + ylab("Component 2") +
    xlim(xl) + ylim(yl)
}




makeTwoPointPlt <- function(X, tmap, lmap, smap, xl, yl, nts = c(0.5, 0.7)) {
# nts <- c(0.5, 0.7)
  nx <- 300
  xgens <- lapply(nts, function(nt) {
    K_y <- lapply(1:2, function(i) cov_matrix(tmap, tmap, as.numeric(lmap[i]), as.numeric(smap[i])))
    K_star <- lapply(1:2, function(i) cov_matrix(tmap, nt, as.numeric(lmap[i])))
    K_dstar <- lapply(1:2, function(i) cov_matrix(nt, nt, as.numeric(lmap[i])))
    
    mu_star <- lapply(1:2, function(i) {
      t(K_star[[i]]) %*% solve(K_y[[i]]) %*% X[,i]
    })
    
    Sigma <- lapply(1:2, function(i) K_dstar[[i]] - t(K_star[[i]]) %*% solve(K_y[[i]]) %*% K_star[[i]] + diag(length(nt)) * smap[i])
    
    Xgen <- sapply(1:2, function(i) rmvnorm(nx, mu_star[[i]], Sigma[[i]]))
    Xgen
  })
  
  
  # construct mu*
  nt <- runif(100)
  K_y <- lapply(1:2, function(i) cov_matrix(tmap, tmap, as.numeric(lmap[i]), as.numeric(smap[i])))
  K_star <- lapply(1:2, function(i) cov_matrix(tmap, nt, as.numeric(lmap[i])))
  K_dstar <- lapply(1:2, function(i) cov_matrix(nt, nt, as.numeric(lmap[i])))
  
  mu_star <- lapply(1:2, function(i) {
    t(K_star[[i]]) %*% solve(K_y[[i]]) %*% X[,i]
  })
  
  
  dd <- data.frame(do.call("rbind", xgens))
  dd$Pseudotime <- rep(c("0.5", "0.7"), each = nrow(dd) / 2)
  
  dx <- data.frame(X)
  dm <- data.frame(do.call("cbind", mu_star))
  dm$t <- nt
  dm <- dplyr::arrange(dm, t)
  
  plt_joint <- ggplot() +
    stat_density2d(data = dd, aes(x = X1, y = X2, fill = Pseudotime, alpha = ..level..), 
                   geom = "polygon", n = 100, contour=TRUE, h = 0.8) +
    geom_point(data = dx, aes(x = X1, y = X2), shape = 21,
               fill = 'black', colour = 'white', size = 3, alpha = 0.6) +
    theme(legend.position="none") + 
    scale_fill_manual (values=c("0.5"="darkred","0.7"="darkblue")) +
    geom_path(data = dm, aes(x = X1, y = X2), size = 1.4, alpha = 0.7) +
    xlab("Component 1") + ylab("Component 2") +
    xlim(xl) + ylim(yl)
  plt_joint
}

set.seed(123)

xlims <- list(c(-2.5, 2), c(-2, 2), c(-2, 2))
ylims <- list(c(-2, 2), c(-4.2, 2.5), c(-2.5, 2.5))

source("gputils/gputils.R")
reps <- c("Xle", "Xle", "Xpca")

post_tracefiles <- paste0("data/", c("trapnell", "burns", "shin"), "_pseudotime_traces.h5")
rep_files <- paste0("data/", c("trapnell", "burns", "shin"), "_embeddings.h5")

## 'diffuseness' plots
plts <- lapply(1:length(post_tracefiles), function(i) {
  post_tracefile <- post_tracefiles[i]
  pst <- h5read(post_tracefile, "pst")
  lambda <- h5read(post_tracefile, "lambda") 
  sigma <- h5read(post_tracefile, "sigma")
  X <- h5read(rep_files[i], reps[i])
  X <- standardize(X)
  
  tmap <- posterior.mode(mcmc(pst))
  lmap <- posterior.mode(mcmc(lambda))
  smap <- posterior.mode(mcmc(sigma))
  makeUncertaintyPlot(X, tmap, lmap, smap, xlims[[i]], ylims[[i]])
})

## two point plots
tplts <- lapply(1:length(post_tracefiles), function(i) {
  post_tracefile <- post_tracefiles[[i]]
  pst <- h5read(post_tracefile, "pst")
  lambda <- h5read(post_tracefile, "lambda") 
  sigma <- h5read(post_tracefile, "sigma")
  X <- h5read(rep_files[i], reps[i])
  X <- standardize(X)
  
  tmap <- posterior.mode(mcmc(pst))
  lmap <- posterior.mode(mcmc(lambda))
  smap <- posterior.mode(mcmc(sigma))
  makeTwoPointPlt(X, tmap, lmap, smap, xlims[[i]], ylims[[i]])
})

all_plots <- c(plts, tplts)
ag <- cowplot::plot_grid(plotlist = all_plots, nrow = 2, labels = c("A", "B", "C", "D", "E", "F"))
outfile <- "figs/envelope/2_cloud.png"
ggsave(ag, filename = outfile, width=8.5, height=6, scale = 1.5)



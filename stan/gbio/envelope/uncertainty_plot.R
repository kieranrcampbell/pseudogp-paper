library(coda)
library(MCMCglmm)
library(mvtnorm)
library(ggplot2)
library(rhdf5)

data(le_fit)
data(monocle_le)

tmcmc <- mcmc(extract(le_fit, pars = "t")$t)
tmap <- posterior.mode(tmcmc)

lmcmc <- mcmc(extract(le_fit, pars = "lambda")[[1]][,1,])
lmap <- posterior.mode(lmcmc)

smcmc <- mcmc(extract(le_fit, pars = "sigma")[[1]][,1,])
smap <- posterior.mode(smcmc)

nnt <- nrow(X)


makeUncertaintyPlot <- function(X, tmap, lmap, smap, nnt = 100, nx = 100, h = 2) {
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
    xlab("Component 1") + ylab("Component 2")
}

base_dir <- "~/mount/GP/pseudogp2/data/"
post_tracefiles <- paste0(base_dir,
                          c("stan_traces_for_gbio.h5", "ear_stan_traces.h5", "waterfall_stan_traces.h5"))

plts <- lapply(post_tracefiles, function(post_tracefile) {
  pst <- h5read(post_tracefile, "pst")
  lambda <- h5read(post_tracefile, "lambda") 
  sigma <- h5read(post_tracefile, "sigma")
  X <- h5read(post_tracefile, "X")
  
  tmap <- posterior.mode(mcmc(pst))
  lmap <- posterior.mode(mcmc(lambda))
  smap <- posterior.mode(mcmc(sigma))
  makeUncertaintyPlot(X, tmap, lmap, smap)
})

pg <- cowplot::plot_grid(plotlist = plts, nrow = 1, labels = c("A", "B", "C"))
ggsave(pg, filename = "2_cloud.png", width = 8.5, height = 3, scale = 1.5)

library(rstan)
library(rhdf5)
library(MCMCglmm)
library(coda)
library(ggplot2)
library(mvtnorm)
library(pscl)
library(dplyr)
library(readr)

plotFromHDF5 <- function(name, h5file) {
  pst <- h5read(h5file, paste0(name, "/pst"))
  lambda <- h5read(h5file, paste0(name, "/lambda"))
  sigma <- h5read(h5file, paste0(name, "/sigma"))
  X <- h5read(h5file, "X")
  return( makeEnvelopePlot(pst, lambda, sigma, X) )
}

## Pinched from envelope.R -------

## plot a bunch of posterior curves
makeEnvelopePlot <- function(pst, lambda, sigma, X, ncurves = 80) {
  set.seed(123)
  ns <- nrow(pst)
  
  pmcs <- lapply(sample(ns, ncurves), function(i) {
    t <- pst[i,]
    l <- lambda[i,]
    s <- sigma[i,]
    posterior_mean_curve(X, t, l, s, nnt = 150)
  })
  
  #mus <- pmc$mu ; nt <- pmc$t
  mus <- lapply(pmcs, function(x) x$mu)
  M <- data.frame(do.call("rbind", mus))
  names(M) <- c("M1", "M2")
  M$curve <- rep(1:ncurves, each = nrow(mus[[1]]))
  M$nt <- unlist(lapply(pmcs, function(x) x$t))
  M <- dplyr::arrange(M, curve, nt)
  
  plt <- ggplot()
  for(i in 1:ncurves) {
    plt <- plt + geom_path(data = filter(M, curve == i), aes(x = M1, y = M2), size = 2, alpha = .07) 
  }
  
  plt <- plt + geom_point(data = data.frame(X), aes(x = X1, y = X2), shape = 21, 
                          fill = 'darkred', colour = 'lightgrey', size = 3, alpha = 0.5) + 
    cowplot::theme_cowplot() + 
    xlab("Component 1") + ylab("Component 2")
  
  return( plt )
}



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



set.seed(123)
base_dir <- "/net/isi-scratch/kieran/"
# base_dir <- "~/mount/"
setwd(file.path(base_dir, "GP/pseudogp2/stan/gbio/varygamma"))
source(file.path(base_dir, "GP/pseudogp2/gputils/gputils.R"))

h5file <- file.path(base_dir, "GP/pseudogp2/data/varygamma_traces.h5")

names <- paste0("g", 1:3)


# posterior curve uncertainty plots ---------------------------------------

eplots <- lapply(names, plotFromHDF5, h5file)
# let's reorder plots in order of shrinkage strength
ep <- eplots[c(1,3,2)]
pg <- cowplot::plot_grid(plotlist = ep, nrow = 1, 
                         labels = c("α = 30, β = 5", "α = 5, β = 1", "α = 10, β = 3"))
ggsave(pg, filename = "S1_varygamma1.png", width = 8, height = 3, scale = 1.5)


# posterior data uncertainty plots ----------------------------------------

plts <- lapply(names, function(name) {
  tmap <- posterior.mode(mcmc(h5read(h5file, paste0(name, "/pst"))))
  lmap <- posterior.mode(mcmc(h5read(h5file, paste0(name, "/lambda"))))
  smap <- posterior.mode(mcmc(h5read(h5file, paste0(name, "/sigma"))))
  X <- h5read(h5file, "X")
  makeUncertaintyPlot(X, tmap, lmap, smap)
})

dplts <- plts[c(1,3,2)]

dg <- cowplot::plot_grid(plotlist = dplts, nrow = 1, 
                         labels = c("α = 30, β = 5", "α = 5, β = 1", "α = 10, β = 3"))
ggsave(dg, filename = "S1_varygamma2.png", width = 8, height = 3, scale = 1.5)





# other stuff -------------------------------------------------------------



gs <- paste0("diffexpr/g", 1:3, ".txt")
g <- lapply(gs, read_csv)

fdrFromDF <- function(df) dplyr::filter(df, test == "ss", type == "falsepos")$value
fdr <- data.frame(GammaAlphaBeta = c("30,5", "10,3", "5,1"),
                  FDR = sapply(g, fdrFromDF))

afdr_plt <- ggplot(fdr) + geom_bar(aes(x = GammaAlphaBeta, y = FDR), stat = "identity") +
  cowplot::theme_cowplot() + 
  theme(legend.position = "none") + coord_flip() +
  ylab("Approximate false discovery rate")

pdf("gamma_fdr.pdf", width = 10, height = 6)
print(pg)
print(afdr_plt)
dev.off()





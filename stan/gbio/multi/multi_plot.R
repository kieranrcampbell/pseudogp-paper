#' Generate envelope plots comparing concurrent fitting to 
#' fitting across different data sources

library(rstan)
library(rhdf5)
library(MCMCglmm)
library(coda)
library(ggplot2)
library(mvtnorm)
library(pscl)


## Pinched from envelope.R -------
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
    plt <- plt + geom_path(data = dplyr::filter(M, curve == i), aes(x = M1, y = M2), size = 2, alpha = .07) 
  }
  
  plt <- plt + geom_point(data = data.frame(X), aes(x = X1, y = X2), shape = 21, 
                          fill = 'darkred', colour = 'lightgrey', size = 3, alpha = 0.5) + 
    cowplot::theme_cowplot() + xlab("") + ylab("") 
  
  return( plt )
}

base_dir <- "~/mount/"
h5file <- paste0(base_dir, "GP/pseudogp2/data/monocle_multi_traces.h5")
h5ls(h5file)


# multi-fit first ---------------------------------------------------------

pst <- h5read(h5file, "multi/pst")

XX <- lapply(c("X", "Y", "Z"), function(x) h5read(h5file, paste0("multi/", x)))
lambda_mat <- h5read(h5file, "multi/lambda") 
lambda <- lapply(list(c(1,2), c(3,4), c(5,6)), function(ind) lambda_mat[,ind])

sigma_mat <- h5read(h5file, "multi/sigma") 
sigma <- lapply(list(c(1,2), c(3,4), c(5,6)), function(ind) sigma_mat[,ind])

multi_plots <- lapply(1:3, function(i) makeEnvelopePlot(pst, lambda[[i]], sigma[[i]], XX[[i]]))
pg <- cowplot::plot_grid(plotlist = multi_plots, nrow = 1, labels = c("Laplacian eigenmaps", "PCA", "t-SNE"))


# individual now ----------------------------------------------------------
h5ls(h5file)
inds <- 1:3
groups <- paste0("indv/g", inds, "/")

extractData <- function(s) lapply(groups, function(g) h5read(h5file, paste0(g, s)))
XXi <- extractData("X")
psti <- extractData("pst")
lambdai <- extractData("lambda")
sigmai <- extractData("sigma")

indv_plots <- lapply(1:3, function(i) makeEnvelopePlot(psti[[i]], lambdai[[i]], sigmai[[i]], XXi[[i]]))
pgi <- cowplot::plot_grid(plotlist = indv_plots, nrow = 1, labels = c("Laplacian eigenmaps", "PCA", "t-SNE"))

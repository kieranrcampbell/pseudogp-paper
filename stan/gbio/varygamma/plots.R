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
    cowplot::theme_cowplot() + xlab("") + ylab("") 
  
  return( plt )
}

set.seed(123)
base_dir <- "/net/isi-scratch/kieran/"
# base_dir <- "~/mount/"
setwd(paste0(base_dir, "GP/pseudogp2/stan/gbio/varygamma"))
source(paste0(base_dir, "GP/pseudogp2/gputils/gputils.R"))

h5file <- paste0(base_dir, "GP/pseudogp2/data/varygamma_traces.h5")

names <- paste0("g", 1:3)

eplots <- lapply(names, plotFromHDF5, h5file)

pg <- cowplot::plot_grid(plotlist = eplots, nrow = 1, 
                         labels = c("α = 30, β = 5", "α = 10, β = 3", "α = 5, β = 1"))
ggsave(pg, filename = "S1_varygamma.png", width = 8, height = 3, scale = 1.5)


gs <- paste0("diffexpr/g", 1:3, ".txt")
g <- lapply(gs, read_csv)

fdrFromDF <- function(df) dplyr::filter(df, test == "ss", type == "falsepos")$value
fdr <- data.frame(GammaAlphaBeta = c("30,5", "10,3", "5,1"),
                  FDR = sapply(g, fdrFromDF))

afdr_plt <- ggplot(fdr) + geom_bar(aes(x = GammaAlphaBeta, y = FDR), stat = "identity") +
  cowplot::theme_cowplot() + 
  theme(legend.position = "none") + coord_flip() +
  ylab("Approximate false discovery rate")



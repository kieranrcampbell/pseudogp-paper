
# This script generates the posterior 'envelope' (uncertainty plots) (unused)

library(ggplot2)
library(reshape2)
library(dplyr)

base_dir <- "~/mount/"

source(file.path(base_dir, "pseudogp-paper/gputils//gputils.R"))

makeEnvelopePlot <- function(post_tracefile) {
  pst <- h5read(post_tracefile, "pst")
  lambda <- h5read(post_tracefile, "lambda") 
  sigma <- h5read(post_tracefile, "sigma")
  X <- h5read(post_tracefile, "X")
  X <- pseudogp::standardize(X)
  
  set.seed(123)
  ns <- nrow(pst)
  
  ncurves <- 80
  
  pmcs <- lapply(sample(ns, ncurves), function(i) {
    t <- pst[i,]
    l <- lambda[i,]
    s <- sigma[i,]
    posterior_mean_curve(X, t, l, s, nnt = 150)
  })
  
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

post_tracefiles <- paste0(base_dir, "pseudogp-paper/",
                    c("monocle_stan_traces.h5", "ear_stan_traces.h5", "waterfall_stan_traces.h5"))

plts <- lapply(post_tracefiles, makeEnvelopePlot)
pg <- cowplot::plot_grid(plotlist = plts, nrow = 1, labels = c("A", "B", "C"))
ggsave(pg, filename = "2_envelope_all.png", width = 8, height = 3, scale = 1.5)
  



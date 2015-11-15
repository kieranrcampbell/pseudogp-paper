
#+ setup, cache=TRUE, message = FALSE
library(rstan)
library(rhdf5)
library(MCMCglmm)
library(coda)
library(ggplot2)
library(mvtnorm)
library(pscl)

set.seed(123)
# base_dir <- "/net/isi-scratch/kieran/"
base_dir <- "~/mount/"
setwd(paste0(base_dir, "GP/pseudogp2/stan/gbio/varygamma"))
source(paste0(base_dir, "GP/pseudogp2/gputils/gputils.R"))

h5file <- paste0(base_dir, "GP/pseudogp2/data/5m_run_with_tau_traces.h5")

X <- h5read(h5file, "X")
X <- apply(X, 2, function(x) (x - mean(x)) / sd(x))
t_gt <- h5read(h5file, "t_gt")

#' We can plot the embedding:
#+ plot-embedding
ggplot(data.frame(X, t_gt)) + geom_point(aes(x = X1, y = X2, color = t_gt), size = 3) +
  scale_color_continuous(low = "darkred", high = "yellow") + theme_bw()

# choose the gamma alpha-betas we're going to use
gab <- data.frame(alpha = c(30, 10, 5), beta = c(5, 3, 1))
gab$mean <- gab$alpha / gab$beta
gab$var <- gab$mean / gab$beta

fits <- list()
for(i in 1:nrow(gab)) {
  data <- list(X = X, N = nrow(X), gamma_alpha = gab$alpha[i], gamma_beta = gab$beta[i])
  fit <- stan(file = "pseudogpvarygamma.stan", data = data, 
               iter = 1000, chains = 1)
  fits <- c(fits, fit)
}

plot(fits[[1]], pars = "t")
plot(fits[[2]], pars = "t")
plot(fits[[3]], pars = "t")

plotFromFit <- function(fit) {
  pst <- extract(fit, "t")$t
  lambda <- extract(fit, "lambda")$lambda
  sigma <- extract(fit, "sigma")$sigma
  makeEnvelopePlot(pst, lambda, sigma, X)
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

eplots <- lapply(fits, plotFromFit)
pg <- cowplot::plot_grid(plotlist = eplots, nrow = 1, 
                         labels = c("α = 30, β = 5", "α = 10, β = 3", "α = 5, β = 1"))
ggsave(pg, filename = "S1_varygamma.png", width = 8, height = 3, scale = 1.5)


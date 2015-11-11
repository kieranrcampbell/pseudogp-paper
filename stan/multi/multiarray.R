library(rstan)
library(rhdf5)
library(MCMCglmm)
library(coda)
library(ggplot2)


base_dir <- "/net/isi-scratch/kieran/"
setwd(paste0(base_dir, "GP/pseudogp2/stan/multi"))
source("../../diffexpr/prep_data.R")
source("../../gputils//gputils.R")


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

# Now individually --------------------------------------------------------

xl <- list(X, Y, Z)
ifits <- lapply(xl, function(X) {
  dx <- array(dim = c(1, 2, nrow(X)))
  dx[1,,] <- t(X)
  data <- list(Ns = 1, P = 2, N = nrow(X), X = dx)
  fit <- stan(file = "pgparray.stan", data = data, 
            iter = 1000, chains = 1)
  return( fit )
})

ifit_posmean_graphs <- lapply(1:3, function(i) {
  fit <- ifits[[i]]
  pst <- extract(fit, "t")
  
  #tmcmc <- mcmc(pst$t)
  tmcmc <- mcmc(pst$t)
  t <- posterior.mode(tmcmc)
  
  lmcmc <- mcmc(extract(fit, "lambda")$lambda[,1,])
  lmap <- posterior.mode(lmcmc)
  
  smcmc <- mcmc(extract(fit, "sigma")$sigma[,1,])
  smap <- posterior.mode(smcmc)

  X <- xl[[i]]
  plot_posterior_mean(X, t_gt, t, lmap, smap)
})

ind_plot <- cowplot::plot_grid(plotlist = ifit_posmean_graphs, ncol = 1, labels = c("Laplacian eigenmaps", "PCA", "t-SNE"))
cowplot::ggsave("ind.png", )

t <- posterior.mode(tmcmc)
colnames(X) <- colnames(Y) <- colnames(X)  <- NULL
xx <- list(X, Y, Z)

plots <- lapply(1:Ns, function(i) plot_posterior_mean(xx[[i]], t_gt, t, lmap[[i]], smap[[i]]))

joint_plot <- plot_grid(plotlist = plots, ncol = 1, labels = c("Laplacian eigenmaps", "PCA", "t-SNE"))
cowplot::ggsave("joint.png", joint_plot, width = 4, height = 5, scale = 2)

all_plot <- plot_grid(joint_plot, ind_plot, ncol = 2, labels = c("Joint", "Individual"))
cowplot::ggsave("comparison.png", all_plot, width = 8, height = 5, scale = 2)

# ind <- sample(nrow(tmcmc), 4)
# plts <- lapply(ind, function(i) plot_posterior_mean(tmcmc[i,], lmcmc[i,], smcmc[i,]))
# cowplot::plot_grid(plotlist = plts)
# 



# 
# tracefile <- "~/mount/GP/pseudogp2/data/stan_traces.h5"
# h5createFile(tracefile)
# h5write(pst$t, tracefile, "pst")


post_tracefile <- "/net/isi-scratch/kieran/GP/pseudogp2/data/envelope_traces.h5"
h5createFile(post_tracefile)
h5write(pst$t, post_tracefile, "pst")
h5write(as.matrix(lmcmc[[1]]), post_tracefile, "lambda")
h5write(as.matrix(smcmc[[2]]), post_tracefile, "sigma")
h5write(X, post_tracefile, "X")
h5write(t_gt, post_tracefile, "t_gt")

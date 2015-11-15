library(rstan)
library(rhdf5)
library(MCMCglmm)
library(coda)
library(ggplot2)
library(scater)

#' Generates traces for multiple datasources (t-SNE, LE, PCA)
#' as well as single and saves to hdf5

#base_dir <- "/net/isi-scratch/kieran/"
base_dir <- "~/mount/"
setwd(paste0(base_dir, "GP/pseudogp2/stan/gbio/multi"))
source("../../diffexpr/monocle/prep_data.R")
source("../../../gputils//gputils.R")

h5outfile <- paste0(base_dir, "GP/pseudogp2/data/monocle_multi_traces.h5")

h5file <- paste0(base_dir, "GP/pseudogp2/data/5m_run_with_tau_traces.h5")


# Laplacian eigenmaps representation --------------------------------------
X <- h5read(h5file, "X") 
X <- apply(X, 2, function(x) (x - mean(x)) / sd(x)) # this is the LE representation
t_gt <- h5read(h5file, "t_gt") # principal curves rep

# PCA representation ------------------------------------------------------
sce <- load_data("~/mount/")
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

if(FALSE) {
fit <- stan(file = "pgparray.stan", data = data, 
            iter = 1000, chains = 1)


pst <- extract(fit, "t")$t

post_map <- posterior.mode(mcmc(pst))

lambda <- lapply(1:Ns, function(i) extract(fit, "lambda")$lambda[,i,])
sigma <- lapply(1:Ns, function(i) extract(fit, "sigma")$sigma[,i,])

lmcmc <- lapply(lambda, mcmc)
smcmc <- lapply(sigma, mcmc)

lmap <- lapply(lmcmc, posterior.mode)
smap <- lapply(lmcmc, posterior.mode)

#plot_posterior_mean(X, t_gt, post_map, lmap[[3]], smap[[3]])

# save to hdf5 ------------------------------------------------------------
if(!file.exists(h5outfile)) h5createFile(h5outfile)
h5createGroup(h5outfile, "multi")
h5write(pst, h5outfile, "multi/pst")
h5write(X, h5outfile, "multi/X")
h5write(Y, h5outfile, "multi/Y")
h5write(Z, h5outfile, "multi/Z")
h5write(do.call("cbind", lambda), h5outfile, "multi/lambda")
h5write(do.call("cbind", sigma), h5outfile, "multi/sigma")
}
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


# save to hdf5 ------------------------------------------------------------

h5createGroup(h5outfile, "indv")
for(i in 1:3) {
  group <- paste0("indv/g", as.character(i))
  h5createGroup(h5outfile, group)
  fit <- ifits[[i]]
  pst <- extract(fit, "t")$t
  
  lambda <- extract(fit, "lambda")$lambda[,1,]
  sigma <- extract(fit, "sigma")$sigma[,1,]

  h5write(xl[[i]], h5outfile, paste0(group, "/X"))
  h5write(pst, h5outfile, paste0(group, "/pst"))
  h5write(lambda, h5outfile, paste0(group, "/lambda"))
  h5write(sigma, h5outfile, paste0(group, "/sigma"))
}

## Done

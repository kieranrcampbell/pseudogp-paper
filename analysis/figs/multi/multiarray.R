library(rstan)
library(rhdf5)
library(MCMCglmm)
library(coda)
library(ggplot2)
library(scater)


# README---- --------------------------------------------------------------
#' KC 15/12/2015 
#' 
#' This should be largely obsolete. Embeddings done in monocle/trajectory_discovery.R
#' and the rest should be re-written as an Rmarkdown doc using library(pseudogp)
#' 
#' Update: working on it 16/12/2015

base_dir <- "~/mount/"
h5outfile <- paste0(base_dir, "GP/pseudogp2/data/monocle_multi_traces.h5")

embedding_files <- paste0(base_dir, "pseudogp-paper/data/",
                          c("trapnell_embeddings.h5", "ear_embeddings.h5", "waterfall_embeddings.h5"))

# for(i in 1:3) {}
i <- 1
ef <- embedding_files[[i]]
reps <- lapply(c("Xle", "Xpca", "Xtsne"), function(r) h5read(ef, r))

smoothing <- list(monocle = c(30, 5), ear = c(9, 1), waterfall = c(8, 2))

indv_fits <- lapply(1:length(reps), function(j) {
  fitPseudotime(reps[[j]], initialise_from = "pca", smoothing_alpha = smoothing[[i]][1], smoothing_beta = smoothing[[i]][2])
})

posteriorCurvePlot(reps[[3]], indv_fits[[3]])

joint_fit <- fitPseudotime(reps, initialise_from = "pca", chains = 2,
                           smoothing_alpha = smoothing[[i]][1], smoothing_beta = smoothing[[i]][2])

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

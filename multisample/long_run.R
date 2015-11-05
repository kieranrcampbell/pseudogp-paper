library(ggplot2)
library(reshape2)
library(rhdf5)
library(coda)
library(MCMCglmm)

h5file <- "/net/isi-scratch/kieran/GP/pseudogp2/data/monocle_gausslik_traces.h5"
h5ls(h5file)

pst <- h5read(h5file, "tchain")
burn_p <- round(nrow(pst)/2)
tmcmc <- mcmc(pst, start = burn_p)

plot(tmcmc)
autocorr.plot(tmcmc, lag.max = 20)

pst_map <- posterior.mode(tmcmc)
t_gt <- h5read(h5file, "t_gt")
plot(t_gt, pst_map)

pst <- data.frame(pst)


pstr <- pst[, sample(1:ncol(pst), 5)]
pstr$iter <- 1:nrow(pstr)

pm <- melt(pstr, id.vars = "iter")

ggplot(pm) + geom_line(aes(x = iter, y = value, color = variable))

lchain <- data.frame(h5read(h5file, "lambda_chain"))
lmcmc <- mcmc(lchain, start = burn_p)
plot(lmcmc)
autocorr.plot(lmcmc, lag.max = 40)

schain <- data.frame(h5read(h5file, "s_chain"))
smcmc <- mcmc(schain, start = burn_p)
plot(smcmc)
autocorr.plot(smcmc, lag.max = 40)



cov_matrix <- function(t1, t2, lambda, sigma = NULL) {
  n1 <- length(t1)
  n2 <- length(t2)
  C <- matrix(NA, nrow = n1, ncol = n2)
  for(i in 1:n1) {
    for(j in 1:n2) {
      C[i, j] <- exp(-lambda * (t1[i] - t2[j])^2)
    }
  }  
  if(!is.null(sigma)) {
    stopifnot(n1 == n2)
    C <- C + sigma * diag(n1)
  }
  return ( C )
}

index <- sample(burn_p:nrow(pst), 1)
t <- pst[index,] ; lambda <- lchain[index,] ; sigma <- schain[index,]

X <- h5read(h5file, "X")
nt <- runif(80)
K_y <- lapply(1:2, function(i) cov_matrix(t, t, as.numeric(lambda[i]), as.numeric(sigma[i])))
K_star <- lapply(1:2, function(i) cov_matrix(t, nt, as.numeric(lambda[i])))
K_dstar <- lapply(1:2, function(i) cov_matrix(nt, nt, as.numeric(lambda[i])))

mu_star <- lapply(1:2, function(i) {
  t(K_star[[i]]) %*% solve(K_y[[i]]) %*% X[,i]
})

mus <- do.call(cbind, mu_star)

ggplot() + geom_point(data = data.frame(X), aes(x = X1, y = X2)) + 
  geom_path(data = data.frame(mus[order(nt),]), aes(x = X1, y = X2), color = "red")


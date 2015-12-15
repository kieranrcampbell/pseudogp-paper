
#+ setup, cache=TRUE, message = FALSE
saveTracesToHDF5 <- function(fit, h5file, name) {
  print(paste("Saving to file", h5file, ", name", name))
  pst <- extract(fit, "t")$t
  lambda <- extract(fit, "lambda")$lambda
  sigma <- extract(fit, "sigma")$sigma
  
  if(!file.exists(h5file)) h5createFile(h5file)
  h5createGroup(h5file, name)
  h5write(pst, h5file, paste0(name, "/pst"))
  h5write(lambda, h5file, paste0(name, "/lambda"))
  h5write(sigma, h5file, paste0(name, "/sigma"))
}

set.seed(123)
# base_dir <- "/net/isi-scratch/kieran/"
base_dir <- "~/mount/"
setwd(paste0(base_dir, "GP/pseudogp2/stan/gbio/varygamma"))
source(paste0(base_dir, "GP/pseudogp2/gputils/gputils.R"))

h5outfile <- paste0(base_dir, "GP/pseudogp2/data/varygamma_traces.h5")

h5file <- paste0(base_dir, "GP/pseudogp2/data/5m_run_with_tau_traces.h5")

X <- h5read(h5file, "X")
X <- apply(X, 2, function(x) (x - mean(x)) / sd(x))
t_gt <- h5read(h5file, "t_gt")

# choose the gamma alpha-betas we're going to use
gab <- data.frame(alpha = c(30, 10, 5), beta = c(5, 3, 1))
gab$mean <- gab$alpha / gab$beta
gab$var <- gab$mean / gab$beta

h5write(X, h5outfile, "X")

fits <- list()
for(i in 1:nrow(gab)) {
  data <- list(X = X, N = nrow(X), gamma_alpha = gab$alpha[i], gamma_beta = gab$beta[i])
  fit <- stan(file = "pseudogpvarygamma.stan", data = data, 
               iter = 1000, chains = 1, seed = 123)
  saveTracesToHDF5(fit, h5outfile, paste0("g", i))
  fits <- c(fits, fit)
}


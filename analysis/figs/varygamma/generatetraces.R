# This creates multiple traces of the monocle dataset for varying gamma

library(pseudogp)
library(rhdf5)

#+ setup, cache=TRUE, message = FALSE
saveTracesToHDF5 <- function(fit, h5file, name) {
  print(paste("Saving to file", h5file, ", name", name))
  pst <- extract(fit, "t")$t
  lambda <- extract(fit, "lambda")$lambda[,1,]
  sigma <- extract(fit, "sigma")$sigma[,1,]
  
  if(!file.exists(h5file)) h5createFile(h5file)
  h5createGroup(h5file, name)
  h5write(pst, h5file, paste0(name, "/pst"))
  h5write(lambda, h5file, paste0(name, "/lambda"))
  h5write(sigma, h5file, paste0(name, "/sigma"))
}

set.seed(123)
# base_dir <- "/net/isi-scratch/kieran/"
base_dir <- "~/mount/"

h5outfile <- file.path(base_dir, "pseudogp-paper/data/monocle_vary_gamma_traces.h5")
h5file <- file.path(base_dir, "pseudogp-paper/data/trapnell_embeddings.h5")

X <- h5read(h5file, "Xle")
# X <- apply(X, 2, function(x) (x - mean(x)) / sd(x))
# t_gt <- h5read(h5file, "t_gt")

# choose the gamma alpha-betas we're going to use
gab <- data.frame(alpha = c(30, 3, 5), beta = c(5, 1, 1))
gab$mean <- gab$alpha / gab$beta
gab$var <- gab$mean / gab$beta

fits <- list()
for(i in 1:nrow(gab)) {
  fit <- fitPseudotime(X, initialise_from = "pca",
                       smoothing_alpha = gab$alpha[i], smoothing_beta = gab$beta[i],
                       seed = 123)
  saveTracesToHDF5(fit, h5outfile, paste0("g", i))
  fits <- c(fits, fit)
}


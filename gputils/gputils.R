## Useful R Gaussian Process functions
## kieran.campbell@sjc.ox.ac.uk

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

posterior_mean_curve <- function(X, t_gt, t, l, s, nnt = 80, reverse = FALSE) {
  if(reverse) t_gt <- 1 - t_gt
  nt <- runif(nnt)
  K_y <- lapply(1:2, function(i) cov_matrix(t, t, as.numeric(l[i]), as.numeric(s[i])))
  K_star <- lapply(1:2, function(i) cov_matrix(t, nt, as.numeric(l[i])))
  K_dstar <- lapply(1:2, function(i) cov_matrix(nt, nt, as.numeric(l[i])))
  
  mu_star <- lapply(1:2, function(i) {
    t(K_star[[i]]) %*% solve(K_y[[i]]) %*% X[,i]
  })
  
  mus <- do.call(cbind, mu_star)
  return(list( mu = mus, t = nt ))
}

plot_posterior_mean <- function(X, t_gt, t, l, s, nnt = 80, reverse = FALSE,
                                curve_color = NULL) {
  colnames(X) <- NULL
  pmc <- posterior_mean_curve(X, t_gt, t, l, s, nnt, reverse)
  mus <- pmc$mu ; nt <- pmc$t
  pdf <- data.frame(mus[order(nt),], nt = nt[order(nt)])
  plt <- ggplot() + 
    geom_point(data = data.frame(X, t_gt), aes(x = X1, y = X2, color = t_gt), size = 3, alpha = 0.5) + 
    theme_bw()
  if(is.null(curve_color)) {
    plt <- plt + geom_path(data = pdf, aes(x = X1, y = X2, color = nt), size = 2, alpha = .8)
  } else {
    plt <- plt + geom_path(data = pdf, aes(x = X1, y = X2), size = 2, alpha = .8, colour = curve_color)
  }
  return( plt )
}


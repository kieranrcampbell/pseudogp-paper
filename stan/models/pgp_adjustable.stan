/* STAN implementation of pseudogp2 

kieran.campbell@sjc.ox.ac.uk

*/

data {
  int<lower = 1> N;
  matrix[N, 2] X;
  real gamma_alpha; // alpha parameter for gamma hyperprior
  real gamma_beta; // as above, beta
  real t_lower; // lower bound for pseudotimes
  real t_upper; // upper bound for pseudotimes
  real tmean; // mean on pseudotime prior
  real tvar; // variance of pseudotime prior
}

transformed data {
  vector[N] X1; 
  vector[N] X2;
  vector[N] mu;

  for(i in 1:N) mu[i] <- 0;

  X1 <- col(X, 1);
  X2 <- col(X, 2);
}

parameters {
  real<lower = 0> lambda[2];
  real<lower = 0> sigma[2];
  real<lower = t_lower, upper = t_upper> t[N];
  real<lower = 0> g;
}

model {
  matrix[N, N] Sigma1;
  matrix[N, N] Sigma2;

  // off-diagonal
  for(i in 1:(N-1)) {
    for(j in (i+1):N) {
      Sigma1[i,j] <- exp(-lambda[1] * pow(t[i] - t[j], 2));
      Sigma2[i,j] <- exp(-lambda[2] * pow(t[i] - t[j], 2));
      Sigma1[j,i] <- Sigma1[i,j];
      Sigma2[j,i] <- Sigma2[i,j];
    }
  }

  // diagonal
  for(k in 1:N) {
    Sigma1[k,k] <- 1 + sigma[1] + 10e-3; // jitter too
    Sigma2[k,k] <- 1 + sigma[2] + 10e-3;
  }

  g ~ gamma(gamma_alpha, gamma_beta);

  for(i in 1:2) {
    sigma[i] ~ inv_gamma(1.0, 1.0);
    lambda[i] ~ exponential(g);
  }

  for(i in 1:N) {
    t[i] ~ normal(tmean, tvar);
  }


  X1 ~ multi_normal(mu, Sigma1);
  X2 ~ multi_normal(mu, Sigma2);
}

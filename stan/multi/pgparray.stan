/* 
STAN implementation of pseudogp2 
This uses multiple sources of data to construct pseudotemporal ordering.

Generalising to any number of data sources of any dimension (but all must be same dimension)

kieran.campbell@sjc.ox.ac.uk
*/

data {
  int<lower = 1> Ns; // number of sources
  int<lower = 1> P; // number of dimensions of each source
  int<lower = 1> N; // 
  vector[N] X[Ns, P]; // Embeddings: an Ns by P array of vector embeddings
}

transformed data {
  vector[N] mu;
  for(i in 1:N) mu[i] <- 0;
}

parameters {
  real<lower = 0> lambda[Ns, P];
  real<lower = 0> sigma[Ns, P];

  real<lower = 0> g[Ns]; // different hyperparameter on lambda for each data source
  real<lower = 0, upper = 1> t[N];
}

model {
  matrix[N, N] Sigma[Ns, P]; // Ns*P covariance matrices

  // off-diagonal
  for(k in 1:Ns) {
    for(l in 1:P) {
      // off-diagonal
      for(i in 1:(N-1)) {
          for(j in (i+1):N) {
            Sigma[k,l,i,j] <- exp(-lambda[k,l] * pow(t[i] - t[j], 2));
            Sigma[k,l,j,i] <- Sigma[k,l,i,j];
        }
      }
      // diagonal
      for(m in 1:N) {
        Sigma[k,l,m,m] <- 1 + sigma[k,l] + 10e-3;
      }
    }
  }

  g ~ gamma(10.0, 3.0);

  for(i in 1:Ns) {
    for(j in 1:P) {
      sigma[i,j] ~ inv_gamma(1.0, 1.0);  
      lambda[i,j] ~ exponential(g[i]);

    }
  }

  for(i in 1:N) {
    t[i] ~ normal(0.5, 1);
  }

  for(i in 1:Ns) {
    for(j in 1:P) {
      X[i,j] ~ multi_normal(mu, Sigma[i,j]);
    }
  }

}

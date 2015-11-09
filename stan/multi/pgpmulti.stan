/* 
STAN implementation of pseudogp2 
This uses multiple sources of data to construct pseudotemporal ordering.
kieran.campbell@sjc.ox.ac.uk

*/

data {
  int<lower = 1> N;
  matrix[N, 2] X; // data type 1
  matrix[N, 2] Y; // data type 2
}

transformed data {
  vector[N] X1; 
  vector[N] X2;
  vector[N] Y1;
  vector[N] Y2;
  vector[N] mu;

  for(i in 1:N) mu[i] <- 0;

  X1 <- col(X, 1);
  X2 <- col(X, 2);
  Y1 <- col(Y, 1);
  Y2 <- col(Y, 2);
}

parameters {
  real<lower = 0> lambdaX[2];
  real<lower = 0> sigmaX[2];
  real<lower = 0> lambdaY[2];
  real<lower = 0> sigmaY[2];
  real<lower = 0, upper = 1> t[N];
  real<lower = 0> g[2];
}

model {
  matrix[N, N] SigmaX1;
  matrix[N, N] SigmaX2;
  matrix[N, N] SigmaY1;
  matrix[N, N] SigmaY2;

  // off-diagonal
  for(i in 1:(N-1)) {
    for(j in (i+1):N) {
      // X cov matrix
      SigmaX1[i,j] <- exp(-lambdaX[1] * pow(t[i] - t[j], 2));
      SigmaX2[i,j] <- exp(-lambdaX[2] * pow(t[i] - t[j], 2));
      SigmaX1[j,i] <- SigmaX1[i,j];
      SigmaX2[j,i] <- SigmaX2[i,j];
  
      // Y cov matrix
      SigmaY1[i,j] <- exp(-lambdaY[1] * pow(t[i] - t[j], 2));
      SigmaY2[i,j] <- exp(-lambdaY[2] * pow(t[i] - t[j], 2));
      SigmaY1[j,i] <- SigmaY1[i,j];
      SigmaY2[j,i] <- SigmaY2[i,j];
    }
  }

  // diagonal
  for(k in 1:N) {
    SigmaX1[k,k] <- 1 + sigmaX[1] + 10e-3; // jitter too
    SigmaX2[k,k] <- 1 + sigmaX[2] + 10e-3;
    SigmaY1[k,k] <- 1 + sigmaY[1] + 10e-3; // jitter too
    SigmaY2[k,k] <- 1 + sigmaY[2] + 10e-3;
  }

  g ~ gamma(10.0, 3.0);

  for(i in 1:2) {
    sigmaX[i] ~ inv_gamma(1.0, 1.0);
    sigmaY[i] ~ inv_gamma(1.0, 1.0);
    lambdaX[i] ~ exponential(g[1]);
    lambdaY[i] ~ exponential(g[2]);
  }

  for(i in 1:N) {
    t[i] ~ normal(0.5, 1);
  }


  X1 ~ multi_normal(mu, SigmaX1);
  X2 ~ multi_normal(mu, SigmaX2);

  Y1 ~ multi_normal(mu, SigmaY1);
  Y2 ~ multi_normal(mu, SigmaY2);
}

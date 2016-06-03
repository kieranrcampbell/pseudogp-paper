## Fitting probabilistic trajectories to Burns et al dataset
### Kieran Campbell <kieran.campbell@sjc.ox.ac.uk>

This document goes through fitting the probabilistic curve to the Burns et al.
dataset. To generate this document in markdown, run `knitr::spin("pseudogp_trapnell.R")` 
from
within R.


```r
library(rstan)
library(rhdf5)
library(MCMCglmm)
library(coda)
library(ggplot2)

devtools::load_all("~/oxford/pseudogp") # TODO - install and replace with `library(pseudogp)`
```

```
## Loading pseudogp
```

```r
set.seed(123)
```

Now read in the data



```r
base_dir <- "~/mount"

h5file <- file.path(base_dir, "pseudogp-paper/data/trapnell_embeddings.h5")
output_hdf5 <- file.path(base_dir, "pseudogp-paper/data/monocle_stan_traces.h5")

X <- h5read(h5file, "Xle")
t_gt <- h5read(h5file, "t_gt")
```

Now fit the probabilistic pseudotime:


```r
fit <- fitPseudotime(X, smoothing_alpha = 30, smoothing_beta = 5, seed = 123)
```

```
## 
## SAMPLING FOR MODEL 'pseudogp' NOW (CHAIN 1).
## 
## Chain 1, Iteration:   1 / 1000 [  0%]  (Warmup)
## Chain 1, Iteration: 100 / 1000 [ 10%]  (Warmup)
## Chain 1, Iteration: 200 / 1000 [ 20%]  (Warmup)
## Chain 1, Iteration: 300 / 1000 [ 30%]  (Warmup)
## Chain 1, Iteration: 400 / 1000 [ 40%]  (Warmup)
## Chain 1, Iteration: 500 / 1000 [ 50%]  (Warmup)
## Chain 1, Iteration: 501 / 1000 [ 50%]  (Sampling)
## Chain 1, Iteration: 600 / 1000 [ 60%]  (Sampling)
## Chain 1, Iteration: 700 / 1000 [ 70%]  (Sampling)
## Chain 1, Iteration: 800 / 1000 [ 80%]  (Sampling)
## Chain 1, Iteration: 900 / 1000 [ 90%]  (Sampling)
## Chain 1, Iteration: 1000 / 1000 [100%]  (Sampling)
## #  Elapsed Time: 147.986 seconds (Warm-up)
## #                72.4939 seconds (Sampling)
## #                220.48 seconds (Total)
```

Plot posterior mean curves


```r
posteriorCurvePlot(X, fit, posterior_mean = TRUE)
```

```
## Plotting traces for 1 representations and 1 chains
```

![plot of chunk posteriorcplt](figure/posteriorcplt-1.png) 

Posterior boxplot of pseudotime distribution


```r
posteriorBoxplot(fit)
```

![plot of chunk posteriorbplt](figure/posteriorbplt-1.png) 

Because we're dealing with stan objects we can very easily extract the posterior samples:


```r
pst <- extract(fit, "t")
tmcmc <- mcmc(pst$t)

smcmc <- mcmc(extract(fit, "sigma")$sigma[,1,]) # need to slice middle index to get single rep.
lmcmc <- mcmc(extract(fit, "lambda")$lambda[,1,])
```

...and save them to HDF5


```r
if(!file.exists(output_hdf5)) h5createFile(output_hdf5)
```

```
## [1] TRUE
```

```r
h5write(X, output_hdf5, "X")
h5write(pst$t, output_hdf5, "pst")
h5write(as.matrix(smcmc), output_hdf5, "sigma")
h5write(as.matrix(lmcmc), output_hdf5, "lambda")
```


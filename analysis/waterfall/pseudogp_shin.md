# Probabilsitic trajectory fitting for the Shin et al. (2015) dataset
### Kieran Campbell <kieran.campbell@sjc.ox.ac.uk>

To turn this into markdown, run `knitr::spin("pseudogp_shin.R")` from
within R.



```r
library(rstan)
library(rhdf5)
library(MCMCglmm)
library(coda)
library(ggplot2)
library(moments)

base_dir <- "~/mount/"

h5file <- file.path(base_dir, "pseudogp-paper/data/waterfall_embeddings.h5")
output_hdf5 <- file.path(base_dir, "pseudogp-paper/data/waterfall_stan_traces.h5")

devtools::load_all("~/oxford/pseudogp") # TODO: change to library(pseudogp)
```

```
## Loading pseudogp
```

```r
set.seed(123)
```

Read in the data


```r
X <- h5read(h5file, "Xpca")
wpst <- h5read(h5file, "wpst") # the pseudotime assigned using Waterfall
```

Fit the pseudotime model


```r
fit <- fitPseudotime(X, initialise_from = "pca", smoothing_alpha = 8, smoothing_beta = 2, seed = 123)
```

```
## Creating pseudotime model with 100 cells and 1 representation(s)
## Standardizing input data
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
## #  Elapsed Time: 37.3406 seconds (Warm-up)
## #                24.591 seconds (Sampling)
## #                61.9316 seconds (Total)
```

Plot the posterior mean curve


```r
posteriorCurvePlot(X, fit)
```

```
## Plotting traces for 1 representations and 1 chains
```

![plot of chunk posmean-curve](figure/posmean-curve-1.png) 

And the boxplots


```r
posteriorBoxplot(fit)
```

![plot of chunk posmean-boxplot](figure/posmean-boxplot-1.png) 

We can extract the posterior traces and compare the MAP estimate to the waterfall fit


```r
pst <- extract(fit, "t")
tmcmc <- mcmc(pst$t)
post_mean <- posterior.mode(tmcmc)
qplot(wpst, post_mean) + theme_bw() + xlab("Waterfall fit") + ylab("Pseudogp map")
```

![plot of chunk compare-map](figure/compare-map-1.png) 

```r
smcmc <- mcmc(extract(fit, "sigma")$sigma[,1,])
lmcmc <- mcmc(extract(fit, "lambda")$lambda[,1,])
```

And finally save everything to HDF5


```r
h5createFile(output_hdf5)
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


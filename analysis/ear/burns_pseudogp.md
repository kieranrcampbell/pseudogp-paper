## Fitting probabilistic trajectories to Burns et al dataset
### Kieran Campbell <kieran.campbell@sjc.ox.ac.uk>

This document goes through fitting the probabilistic curve to the Burns et al.
dataset. To generate this document in markdown, run `knitr::spin("burns_pseudogp.R")` 
from
within R.
Quick setup



```r
library(rstan)
library(rhdf5)
library(MCMCglmm)
library(coda)
library(ggplot2)
library(moments)
library(devtools)

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

h5file <- file.path(base_dir, "pseudogp-paper/data/ear_embeddings.h5")
output_hdf5 <- file.path(base_dir, "pseudogp-paper/data/ear_stan_traces.h5")

X <- h5read(h5file, "Xle")
X <- apply(X, 2, function(x) (x - mean(x)) / sd(x))
t_gt <- h5read(h5file, "t_gt")
```

and quickly plot to make sure it looks right


```r
ggplot(data.frame(X, t_gt)) + 
  geom_point(aes(x = X1, y = X2, color = t_gt)) + theme_bw()
```

![plot of chunk quick-plot](figure/quick-plot-1.png) 

### Fit the pseudotime


```r
fit <- fitPseudotime(X, initialise_from = "principal_curve", 
                     smoothing_alpha = 9, smoothing_beta = 1, seed = 123)
```

```
## Creating pseudotime model with 95 cells and 1 representation(s)
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
## #  Elapsed Time: 56.1843 seconds (Warm-up)
## #                39.1821 seconds (Sampling)
## #                95.3664 seconds (Total)
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


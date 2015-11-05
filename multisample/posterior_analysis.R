
#' # Convergence analysis of MCMC chains
#' This document looks at (rough measures of) the convergence of MCMC chains
#' of the student-likelihood GP model using 50% subsamples of the data.

#+ setup
library(rhdf5)
library(reshape2)
library(ggplot2)
library(dplyr)
library(modeest)

get_map <- function(x) mlv(x,method = "HSM")$M

data_dir <- "/net/isi-project/CW010_CAMPBELL_SCNGEN/data/GP/multisample"
h5file <- paste0(data_dir, "/multitrace.h5")
nsamples <- 50

#' ### Slope of posterior trace
#' If the posteriors are converged, the trace should be approximately 1

#+ posterior-slopes
slopes <- sapply(1:nsamples, function(i) {
  sname <- paste0("sample_", i, "/tchain") 
  trace <- h5read(h5file, sname)
  burn_period <- as.integer(nrow(trace) / 2)
  approx_slopes <- apply(trace[burn_period:nrow(trace),], 2,
                         function(tr) {
                           iter <- (1:length(tr)) / length(tr) # normalise 0-1 so slopes make sense
                           fit <- lm(tr ~ iter)
                           coef(fit)[2]
                         })
})

#sm <- melt(data.frame(slopes), value.name = "y")
#ggplot(sm) + geom_boxplot(aes(y=y, x=factor(0))) + facet_wrap(~variable)  

#' We can look at boxplots of the slopes:

#+ boxplots
ggplot(data.frame(x = as.vector(slopes))) + geom_boxplot(aes(y=x, x = factor(0)))

#' Clearly there are quite a few outliers. Let's look at a few:

#+ bad-outliers
bad <- which(abs(slopes[,2]) > 0.1)
bpst <- h5read(h5file, "sample_2/tchain")
burn_period <- as.integer(nrow(bpst) / 2)
bpst <- bpst[burn_period:nrow(bpst), bad]
bpst <- data.frame(bpst)
bpst$iter <- 1:nrow(bpst)
mpst <- melt(bpst, value.name = "pseudotime", id.vars = "iter")
ggplot(mpst) + geom_line(aes(x = iter, y = pseudotime)) + facet_wrap(~variable)

#' It appears most of them can be explained by multi-modality.

#' ### Comparison to principal curve fit
#' Load in data and calculate correlations

#+ load-calculate, warning = FALSE
h5data <- "/net/isi-scratch/kieran/GP/pseudogp2/data/embeddings.h5"
pst <- h5read(h5data, "monocle/pseudotime")

cors <- sapply(1:nsamples, function(i) {
  sname <- paste0("sample_", i)
  trace <- h5read(h5file, paste0(sname, "/tchain"))
  to_sample <- h5read("/net/isi-scratch/kieran/GP/pseudogp2/data/multisample.h5", 
                      paste0(sname, "/to_sample"))
  map_pseudotime <- apply(trace[burn_period:nrow(trace),], 2, get_map)
  cor(pst[to_sample], map_pseudotime)
})

ggplot(data.frame(x = 1:length(cors), y = cors)) + geom_point(aes(x = x, y = y), size = 3) +
  theme_bw() + scale_y_continuous(breaks = seq(-.9, .9, length.out = 20)) +
  xlab("Subsample") + ylab("Correlation to posterior mean pseudotime")

#' Now we want to plot different subsamples against the pc fits
#' to identify major sources of variablility

#+ plot-against-pc
plot_against_pc <- function(i) {
  sname <- paste0("sample_", i)
  trace <- h5read(h5file, paste0(sname, "/tchain"))
  to_sample <- h5read("/net/isi-scratch/kieran/GP/pseudogp2/data/multisample.h5", 
                      paste0(sname, "/to_sample"))
  map_pseudotime <- apply(trace[burn_period:nrow(trace),], 2, get_map)
  qplot(pst[to_sample], map_pseudotime)
}

#' Let's look at the really bad points
#+ bad, warning = FALSE
trouble_points <- abs(cors) < 0.6
trouble_plots <- lapply(which(trouble_points), plot_against_pc)
cowplot::plot_grid(plotlist = plot_list)

#' The quite bad ones:
#+ quite-bad, warning = FALSE
boundary_points <- abs(cors) > .6 & abs(cors) < .7
boundary_plots <- lapply(which(boundary_points), plot_against_pc)
cowplot::plot_grid(plotlist = boundary_plots)

#' The okay ones:
#+ okay, warning = FALSE
ok_points <- abs(cors) > .7 & abs(cors) < .8
ok_plots <- lapply(which(ok_points), plot_against_pc)
cowplot::plot_grid(plotlist = ok_plots)

#' And the 'good' ones:
#+ good, warning = FALSE
good_points <- abs(cors) > .8
good_plots <- lapply(which(good_points), plot_against_pc)
cowplot::plot_grid(plotlist = good_plots)

#' Generate this document with
#+ generate-doc
# knitr::spin("posterior_analysis.R")

#' Sessioninfo:
#+ sessioninfo
sessionInfo()


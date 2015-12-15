library(rstan)
library(rhdf5)
library(coda)
library(MCMCglmm)
library(RColorBrewer)



pcurvePlot <- function(pst, lambda, sigma, X) {
  ns <- nrow(pst)

  ncurves <- ncol(pst)
  getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
  cols <- getPalette(ncurves)

  pmcs <- lapply(1:ncurves, function(i) {
    t <- pst[,i]
    l <- lambda[,i]
    s <- sigma[,i]
    posterior_mean_curve(X, t, l, s, nnt = 150)
  })

  #mus <- pmc$mu ; nt <- pmc$t
  mus <- lapply(pmcs, function(x) x$mu)
  M <- data.frame(do.call("rbind", mus))
  names(M) <- c("M1", "M2")
  M$curve <- rep(1:ncurves, each = nrow(mus[[1]]))
  M$nt <- unlist(lapply(pmcs, function(x) x$t))
  M <- dplyr::arrange(M, curve, nt)

  plt <- ggplot()
  for(i in 1:ncurves) {
    plt <- plt + geom_path(data = dplyr::filter(M, curve == i), aes(x = M1, y = M2),
                           size = 2, alpha = 0.3, color = 'black')
  }

  plt <- plt + geom_point(data = data.frame(X), aes(x = X1, y = X2), shape = 21,
                          fill = 'darkred', colour = 'lightgrey', size = 3, alpha = 0.5) +
    cowplot::theme_cowplot() + xlab("") + ylab("")

  return( plt )
}

set.seed(123)
base_dir <- "/net/isi-scratch/kieran/"
#base_dir <- "~/mount/"
setwd(paste0(base_dir, "GP/pseudogp2/stan/gbio/chains"))
source(paste0(base_dir, "GP/pseudogp2/gputils/gputils.R"))

h5outfile <- paste0(base_dir, "GP/pseudogp2/data/chains.h5")

h5file <- paste0(base_dir, "GP/pseudogp2/data/5m_run_with_tau_traces.h5")

X <- h5read(h5file, "X")
X <- apply(X, 2, function(x) (x - mean(x)) / sd(x))
t_gt <- h5read(h5file, "t_gt")

# choose the gamma alpha-betas we're going to use
gab <- data.frame(alpha = c(30, 3), beta = c(5, 1))
gab$mean <- gab$alpha / gab$beta
gab$var <- gab$mean / gab$beta

h5write(X, h5outfile, "X")

fits <- list()

for(i in 1:nrow(gab)) {
  data <- list(X = X, N = nrow(X), gamma_alpha = gab$alpha[i], gamma_beta = gab$beta[i])
  fit <- stan(file = "pseudogpvarygamma.stan", data = data,
               iter = 1000, chains = 20, seed = 123)
  fits <- c(fits, fit)
}

save(fits, file = "fits.Rdata", compress = "gzip")

fit1 <- fits[[1]]
fit2 <- fits[[2]]

gen_plots <- lapply(fits, function(fit) {
  pst <- extract(fit, pars = "t", permute = FALSE)
  lambda1 <- extract(fit, pars = "lambda", permute = FALSE)
  sigma1 <- extract(fit, pars = "sigma", permute = FALSE)
  t_modes <- apply(pst, 2, function(x) posterior.mode(mcmc(x)))
  lm <- apply(lambda1, 2, function(x) posterior.mode(mcmc(x)))
  sm <- apply(sigma1, 2, function(x) posterior.mode(mcmc(x)))


  t_hpd <- apply(pst, 2, function(x) {
    hpd <- HPDinterval(mcmc(x), prob = 0.95)
    return(hpd[,2] - hpd[,1])
    })

  t_mode_cor <- cor(t_modes)
  lt <- lower.tri(t_mode_cor)
  lt[lt == FALSE] <- NA
  ltc <- as.vector(lt * abs(t_mode_cor))
  ltc <- ltc[!is.na(ltc)]

  hpd_cor <- cor(t_hpd)
  hc <- as.vector(lt * abs(hpd_cor))
  hc <- hc[!is.na(hc)]

  pcplt <- pcurvePlot(t_modes, lm, sm, X)

  mapplt <- qplot(ltc) + cowplot::theme_cowplot() +
    xlab("Correlation between MAP pseudotime estimates")

  hpdplt <- qplot(hc) + cowplot::theme_cowplot() +
    xlab("Correlation between HPD interval widths")

  return(list(pcplt = pcplt, mapplt = mapplt, hpdplt = hpdplt))
})

pdf("vary_lambda_chain.pdf", width = 12, height = 7)
cowplot::plot_grid(gen_plots[[1]]$pcplt, gen_plots[[2]]$pcplt,
                   labels = c("mean = 6, var = 1.2", "mean = 3, var = 3"))
cowplot::plot_grid(gen_plots[[1]]$mapplt, gen_plots[[2]]$mapplt,
                   labels = c("mean = 6, var = 1.2", "mean = 3, var = 3"))
cowplot::plot_grid(gen_plots[[1]]$hpdplt, gen_plots[[2]]$hpdplt,
                   labels = c("mean = 6, var = 1.2", "mean = 3, var = 3"))
dev.off()


## plots of marker genes

source(paste0(base_dir, "GP/pseudogp2/stan/diffexpr/monocle/prep_data.R"))
sce <- load_data(base_dir)

cdk1 <- match("CDK1", fData(sce)$gene_short_name)
myog <- match("MYOG", fData(sce)$gene_short_name)

redDim(sce) <- X
cdk1_plt <- plotReducedDim(sce, size_by = featureNames(sce)[cdk1])
myog_plt <- plotReducedDim(sce, size_by = featureNames(sce)[myog])

gplt <- cowplot::plot_grid(cdk1_plt, myog_plt, nrow = 1, labels = c("CDK1", "MYOG"))

ggsave("S1b_marker_genes.png", gplt, width = 8, height = 3, scale = 1.5)

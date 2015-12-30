library(rstan)
library(rhdf5)
library(coda)
library(MCMCglmm)
library(RColorBrewer)
library(pseudogp)
library(cowplot)

set.seed(123)
base_dir <- "/net/isi-scratch/kieran/"
#base_dir <- "~/mount/"

#setwd(paste0(base_dir, "GP/pseudogp2/stan/gbio/chains"))
#source(paste0(base_dir, "GP/pseudogp2/gputils/gputils.R"))

h5outfile <- file.path(base_dir, "pseudogp-paper/data/chains.h5")

embedding_file <- file.path(base_dir, "pseudogp-paper/data/trapnell_embeddings.h5")

X <- h5read(embedding_file, "Xle")
t_gt <- h5read(embedding_file, "t_gt")

# choose the gamma alpha-betas we're going to use
gab <- data.frame(alpha = c(35, 3), beta = c(5, 1))
gab$mean <- gab$alpha / gab$beta
gab$var <- gab$mean / gab$beta
nchains <- 10

fits <- list()

for(i in 1:nrow(gab)) {
  gamma_alpha = gab$alpha[i]
  gamma_beta = gab$beta[i]
  fit <- fitPseudotime(X, smoothing_alpha = gamma_alpha, smoothing_beta = gamma_beta,
                       chains = nchains)
  fits <- c(fits, fit)
}

save(fits, file = file.path(base_dir, "pseudogp-paper/analysis/figs/chains/fits.Rdata"), compress = "gzip")
load(file.path(base_dir, "pseudogp-paper/analysis/figs/chains/fits.Rdata"))

pcplts <- lapply(1:2, function(i) posteriorCurvePlot(X, fits[[i]]))
# plot_grid(plotlist = pcplts, nrow = 1, labels = c("A", "B"))

celln <- 100
x <- sapply(1:2, function(i) extract(fits[[i]], "t")$t)
hplts <- lapply(1:2, function(i) ggplot(data.frame(t = x[,i])) + geom_density(aes(x = t), fill = "black"))

gplt <- plot_grid(plotlist = c(pcplts, hplts), labels = c("A", "B", "C", "D"))

ggsave(file.path(base_dir, "pseudogp-paper/analysis/figs/chains/chains.png"), gplt,
       width = 8, height = 5)



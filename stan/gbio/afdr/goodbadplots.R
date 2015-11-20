
library(dplyr)
library(ggplot2)
library(readr)
library(RColorBrewer)
library(grid)
library(rhdf5)
library(embeddr)
library(coda)
library(MCMCglmm)

#base_dir <- "/net/isi-scratch/kieran/"

makeGXPlot <- function(gene, sce, pst, to_fit = 200) {

  inds <- sample(nrow(pst), to_fit)
  
  gi <- grep(gene, fData(sce)$gene_short_name) 
  
  models <- apply(pst[inds,], 1, function(t) {
    sce$pseudotime <- t
    fit <- fit_pseudotime_model(sce, gi) 
    return( predict(fit) )
  })
  
  
  models[models < sce@lowerDetectionLimit] <- sce@lowerDetectionLimit
  
  post_mean <- posterior.mode(mcmc(pst))
  
  sce$pseudotime <- post_mean
  map_model <- predict(fit_pseudotime_model(sce, gi))
  map_model[map_model < sce@lowerDetectionLimit] <- sce@lowerDetectionLimit
  
  y <- exprs(sce)[gi,]
  
  dfx <- data.frame(t = post_mean, y = y)
  
  plt <- ggplot() + geom_point(data = dfx, aes(x = t, y = y), shape = 21, 
                               fill = 'darkred', colour = 'lightgrey', size = 3, alpha = 0.5)
  
  for(i in 1:to_fit) {
    plt <- plt + geom_line(aes(x = x, y = y), data = data.frame(x = pst[inds[i],], y = models[,i]),
                  size = 2, alpha = 0.03)
  }
  
  plt <- plt + xlab("Pseudotime") + ylab("Expression") +
    cowplot::theme_cowplot()
   plt <- plt + geom_line(aes(x = x, y = y), data = data.frame(x = post_mean, y = map_model),
                         size = 1, color = "red")
  
  return(plt)
}

base_dir <- "~/mount/"
setwd(paste0(base_dir, "GP/pseudogp2/stan/gbio/afdr"))
source(paste0(base_dir, "GP/pseudogp2/stan/diffexpr/monocle/prep_data.R"))
pstfile <- paste0(base_dir, "GP/pseudogp2/data/monocle_multi_traces.h5")

sce <- load_data(base_dir)
pst <- h5read(pstfile, "multi/pst")



set.seed(123)
genes <- c("ITGAE", "^ID1$")
plts <- lapply(genes, makeGXPlot, sce, pst)

gridplt <- cowplot::plot_grid(plotlist = plts, nrow = 1, labels = c("ITGAE", "ID1"))
ggsave("goodbadplots.png", gridplt, width = 4, height = 2, scale = 1.3)





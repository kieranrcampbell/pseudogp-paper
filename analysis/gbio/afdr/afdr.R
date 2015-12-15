
library(dplyr)
library(ggplot2)
library(readr)
library(RColorBrewer)
library(grid)
library(rhdf5)
library(embeddr)
library(coda)

base_dir <- "/net/isi-scratch/kieran/"
setwd(paste0(base_dir, "GP/pseudogp2/stan/gbio/workflow"))
source(paste0(base_dir, "GP/pseudogp2/stan/diffexpr/monocle/prep_data.R"))
pstfile <- paste0(base_dir, "GP/pseudogp2/data/stan_traces_for_gbio.h5")

sce <- load_data()
pst <- h5read(pstfile, "pst")

to_fit <- 200
set.seed(123)
inds <- sample(nrow(pst), to_fit)

gene <- "^ID1$"
gi <- grep(gene, fData(sce)$gene_short_name) 

models <- apply(pst[inds,], 1, function(t) {
  sce$pseudotime <- t
  fit <- fit_pseudotime_model(sce, gi) 
  return( predict(fit) )
})

models[models < sce@lowerDetectionLimit] <- sce@lowerDetectionLimit

post_mean <- posterior.mode(mcmc(pst))

y <- exprs(sce)[gi,]

dfx <- data.frame(t = post_mean, y = y)

plt <- ggplot() + geom_point(data = dfx, aes(x = t, y = y), shape = 21, 
                             fill = 'darkred', colour = 'lightgrey', size = 3, alpha = 0.5)

for(i in 1:to_fit) {
  plt <- plt + geom_line(aes(x = x, y = y), data = data.frame(x = pst[inds[i],], y = models[,i]),
                size = 2, alpha = 0.07)
}

plt <- plt + xlab("Pseudotime") + ylab("Expression")

ggsave("multiexpr.png", plt, width = 4, height = 3, scale = 1.3)





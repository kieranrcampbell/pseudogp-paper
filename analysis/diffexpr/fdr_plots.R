## Generate all FDR plots
## Edit 18/12/2015 - this is now all figure 5

library(dplyr)
library(ggplot2)
library(readr)
library(RColorBrewer)
library(grid)
library(MCMCglmm)
library(coda)
library(cowplot)
library(gplots)
library(grid)

base_dir <- "~/mount"
fdr_dir <- file.path(base_dir, "pseudogp-paper/analysis/diffexpr")

source(file.path(fdr_dir, "common.R"))

monocle <- read_csv(file.path(fdr_dir, "monocle_fdr.txt"))
ear <- read_csv(file.path(fdr_dir, "ear_fdr.txt"))
waterfall <- read_csv(file.path(fdr_dir, "waterfall_fdr.txt"))

fdrFromDF <- function(df) filter(df, test == "ss", type == "falsepos")$value
fdr <- data.frame(Publication = c("Trapnell 2014", "Burns 2015", "Shin 2015"),
                  FDR = sapply(list(monocle, ear, waterfall), fdrFromDF))


afdr_plt <- ggplot(fdr) + geom_bar(aes(x = Publication, y = FDR), stat = "identity") +
  cowplot::theme_cowplot() + 
  theme(legend.position = "none") + coord_flip() +
  ylab("Approximate false discovery rate")
# scale_fill_manual(name = "", values = brewer.pal(3, "Set1") ) +

#ggsave(file.path(fdr_dir, "afdr.png"), afdr_plt, width = 6, height = 5, scale = 1.6)

h5_diffexpr <- file.path(base_dir, "pseudogp-paper/data/monocle_diffexpr.h5")
pstfile <- file.path(base_dir, "pseudogp-paper/data/monocle_stan_traces.h5")
sce_file <- file.path(base_dir, "pseudogp-paper/data/sce_monocle.Rdata")
load(sce_file)
sce <- sce_23

n_cells_exprs <- rowSums(exprs(sce) > sce@lowerDetectionLimit)
genes_to_use <- n_cells_exprs > (0.1 * ncol(sce)) # select genes expressed in at least 10% of cells
sce <- sce[genes_to_use,]

pst <- h5read(pstfile, "pst")

sigList <- pvalsFromHDF5(h5_diffexpr)

pst_map <- posterior.mode(mcmc(pst))
sce$pseudotime <- pst_map

de_test <- pseudotime_test(sce, n_cores = 1)

dfc <- data.frame(psig = sigList$ss_prop_sig, qval = de_test$q_val)
dfc <- mutate(dfc, is_sig = qval < 0.05)

assignT <- function(r) {
  if(r[1] >= 0.95 && r[2] < 0.05) {
    return("TP")
  } else if(r[1] < 0.95 && r[2] < 0.05) {
    return("FP")
  } else if(r[1] < 0.95 && r[2] >= 0.05) {
    return("TN") 
  } else if(r[1] >= 0.95 && r[2] >= 0.05) {
    return("FN")
  }
}

dfc$class <- apply(dfc, 1, assignT)

fe_plt <- ggplot(dfc) + geom_point(aes(x = psig, y = qval, color = class), alpha = 0.5, size = 2.5) +
  scale_color_manual(name = "", labels = c("False Positive", "True Negative", "True Positive"),
                    values = brewer.pal(4, "Set1")) +
  xlab("Proportion significant") + ylab("FDR-adjusted Q-value") +
  geom_vline(xintercept = 0.95, linetype = 2, alpha = 0.5) +
  geom_hline(yintercept = 0.05, linetype = 2, alpha = 0.5) +
  theme(legend.position = c(0.7,0.7))

# ggsave("fe.png", fe_plt, width = 5, height = 3, scale = 1.5)

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
  plt <- plt + scale_x_continuous(breaks = c(0, 0.5, 1))
  
  return(plt)
}



set.seed(123)
genes <- c("ITGAE", "^ID1$")
plts <- lapply(genes, makeGXPlot, sce, pst)

gridplt <- cowplot::plot_grid(plotlist = plts, nrow = 1, labels = c("ITGAE", "ID1"), 
                              hjust = c(-2, -4), vjust = c(0, 0), label_size = 12)
#ggsave("goodbadplots.png", gridplt, width = 4, height = 2, scale = 1.3)


# heatmap -----------------------------------------------------------------

makeHeatmapPlot <- function(gene, sce, pst) {
  set.seed(123)
  gi <- grep(gene, fData(sce)$gene_short_name) 
  x <- exprs(sce[gi,])
  to_sample <- 20
  psts_sampled <- pst[sample(1:nrow(pst), to_sample, replace=FALSE),]
  orders <- apply(psts_sampled, 1, order)
  gex <- apply(orders, 2, function(o) x[o])
  
  df <- data.frame(gex)
  df$pseudotime_order <- 1:nrow(df)
  dfm <- melt(df, id.vars = "pseudotime_order")
  ggplot(dfm) + geom_tile(aes(x = pseudotime_order, y = variable, fill = value)) +
    scale_fill_gradient(low = "white", high = "steelblue") +
    theme(legend.position = "none", axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
    theme(axis.line = element_blank()) +
    xlab("Pseudotime order") + ylab(expression("Posterior\nsample")) 
}

heatplots <- lapply(genes, makeHeatmapPlot, sce, pst)
heatplt <- plot_grid(plotlist = heatplots, nrow = 2, 
                     labels = c("ITGAE", "ID1"), label_size = 12, hjust = c(-2, -4)) 

plots <- list(fe_plt, afdr_plt, gridplt, heatplt)
# add margins
plots <- lapply(plots, function(plt) plt + theme(plot.margin = unit(rep(1, 4), "cm")))

# lower_grid <- plot_grid(afdr_plt, gridplt, nrow = 1, labels = c("B", "C"), label_size = 16)

plots <- plots[c(3,4,1,2)]
all_plt <- cowplot::plot_grid(plotlist = plots, nrow = 2, 
                              labels = c("A", "B", "C", "D"), label_size = 16, rel_heights = c(2.5,3))

ggsave(file.path(fdr_dir, "fdr.png"), all_plt, width = 8, height = 6, scale = 1.5)



## Generate all FDR plots

library(dplyr)
library(ggplot2)
library(readr)
library(RColorBrewer)
library(grid)

source("common.R")
source(paste0(base_dir, "GP/pseudogp2/stan/diffexpr/monocle/prep_data.R"))

base_dir <- "/net/isi-scratch/kieran/"
setwd(paste0(base_dir, "GP/pseudogp2/stan/gbio/multi/diffexpr"))

joint <- read_csv("multi_fdr.txt")
le <- read_csv("le_fdr.txt")
pca <- read_csv("pca_fdr.txt")
tsne <- read_csv("tsne_fdr.txt")

fdrFromDF <- function(df) dplyr::filter(df, test == "ss", type == "falsepos")$value
fdr <- data.frame(Representation = c("Joint", "Laplacian eigenmaps", "PCA", "t-SNE"),
                  FDR = sapply(list(joint, le, pca, tsne), fdrFromDF))


afdr_plt <- ggplot(fdr) + geom_bar(aes(x = Representation, y = FDR), stat = "identity") +
  cowplot::theme_cowplot() + 
  theme(legend.position = "none") + coord_flip() +
  ylab("Approximate false discovery rate")
# scale_fill_manual(name = "", values = brewer.pal(3, "Set1") ) +

ggsave("rep_afdr.png", afdr_plt, width = 6, height = 5, scale = 1.3)
# 
# h5file <- paste0(base_dir, "GP/pseudogp2/data/monocle_diffexpr.h5")
# pstfile <- paste0(base_dir, "GP/pseudogp2/data/stan_traces_for_gbio.h5")
# pst <- h5read(pstfile, "pst")
# 
# 
# sce <- load_data()
# sigList <- pvalsFromHDF5(h5file)
# 
# pst_map <- posterior.mode(mcmc(pst))
# sce$pseudotime <- pst_map
# 
# de_test <- pseudotime_test(sce, n_cores = 1)
# 
# dfc <- data.frame(psig = sigList$ss_prop_sig, qval = de_test$q_val)
# dfc <- mutate(dfc, is_sig = qval < 0.05)
# 
# assignT <- function(r) {
#   if(r[1] >= 0.95 && r[2] < 0.05) {
#     return("TP")
#   } else if(r[1] < 0.95 && r[2] < 0.05) {
#     return("FP")
#   } else if(r[1] < 0.95 && r[2] >= 0.05) {
#     return("TN") 
#   } else if(r[1] >= 0.95 && r[2] >= 0.05) {
#     return("FN")
#   }
# }
# 
# dfc$class <- apply(dfc, 1, assignT)
# 
# fe_plt <- ggplot(dfc) + geom_point(aes(x = psig, y = qval, color = class), alpha = 0.5, size = 2.5) +
#   scale_color_manual(name = "", labels = c("False Positive", "True Negative", "True Positive"),
#                     values = brewer.pal(4, "Set1")) +
#   xlab("Proportion significant") + ylab("FDR-adjusted Q-value") +
#   geom_vline(xintercept = 0.95, linetype = 2, alpha = 0.5) +
#   geom_hline(yintercept = 0.05, linetype = 2, alpha = 0.5) +
#   theme(legend.position = c(0.7,0.7))
# 
# ggsave("fe.png", fe_plt, width = 5, height = 3, scale = 1.5)
# 



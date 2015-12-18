library(rstan)
library(rhdf5)
library(MCMCglmm)
library(coda)
library(ggplot2)
library(scater)
library(cowplot)

# README---- --------------------------------------------------------------
#' KC 15/12/2015 
#' 
#' This should be largely obsolete. Embeddings done in monocle/trajectory_discovery.R
#' and the rest should be re-written as an Rmarkdown doc using library(pseudogp)
#' 
#' Update: working on it 16/12/2015
#' 
#' Update: I just want it to be christmas already 17/12/2015

plot_one_dataset <- function(i) {
  set.seed(123)
  embedding_files <- paste0(base_dir, "pseudogp-paper/data/",
                            c("trapnell_embeddings.h5", "ear_embeddings.h5", "waterfall_embeddings.h5"))
  ef <- embedding_files[[i]]
  reps <- lapply(c("Xle", "Xpca", "Xtsne"), function(r) h5read(ef, r))
  reps <- lapply(reps, standardize)
  
  smoothing <- list(monocle = c(30, 5), ear = c(18, 2), waterfall = c(8, 2))
  
  indv_fits <- lapply(1:length(reps), function(j) {
    fitPseudotime(reps[[j]], initialise_from = "pca", smoothing_alpha = smoothing[[i]][1], 
                  smoothing_beta = smoothing[[i]][2], seed = 123)
  })
  
  indv_plts <- lapply(1:3, function(i) {
    posteriorCurvePlot(reps[[i]], indv_fits[[i]], posterior_mean = FALSE, curve_alpha = 0.1, 
                       standardize_ranges = TRUE)
  })
  
  indv_grid <- plot_grid(plotlist = indv_plts, nrow = 1)
  
  joint_fit <- fitPseudotime(reps, initialise_from = "pca",
                             smoothing_alpha = smoothing[[i]][1], 
                             smoothing_beta = smoothing[[i]][2],
                             seed = 123)
  
  joint_grid <- posteriorCurvePlot(reps, joint_fit, curve_alpha = 0.1, 
                                   posterior_mean = FALSE, grid_nrow = 1,
                                   standardize_ranges = TRUE)
  
  return(list(indv = indv_grid, joint = joint_grid))
}

base_dir <- "~/mount/"
h5outfile <- paste0(base_dir, "GP/pseudogp2/data/monocle_multi_traces.h5")

all_dataset_plots <- lapply(1:3, plot_one_dataset)
all_dataset_plots[[2]] <- plot_one_dataset(2)



by_dataset <- lapply(all_dataset_plots, function(dp) {
  plot_grid(dp[[1]], dp[[2]], nrow = 2, labels = c("Individual", "Joint"), label_size = 16, vjust = c(1.2,1.2), hjust = c(-0.5, -1))
})

file_names <- c("a_waterfall.png", "b_burns.png", "c_shin.png")
file_paths <- sapply(1:3, function(i) file.path(base_dir, "pseudogp-paper/analysis/figs/multi", file_names[i]))

for(i in 1:3) ggsave(by_dataset[[i]], filename = file_paths[i], width=6, height=3, scale = 1.8)

# dataset_names <- c("Trapnell", "Burns", "Shin")
# plot_grid(plotlist = by_dataset, nrow = 3, labels = dataset_names, 
#           label_size = 16)

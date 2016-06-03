library(ggplot2)
library(reshape2)
library(dplyr)
library(RColorBrewer)
library(grid)
library(coda)
library(matrixStats)
library(pseudogp)
library(rhdf5)


post_tracefile <- "data/trapnell_pseudotime_traces.h5"

pst <- h5read(post_tracefile, "pst")

ind <- c(21, 121, 131, 151)

d <- data.frame(pst[,ind])
names(d) <- c("Cell1", "Cell2", "Cell3", "Cell4")
dm <- melt(d)

cols <- scale_colour_brewer(palette = "Set1", type = "qual")

plt <- ggplot(dm) + geom_density(aes(x = value, fill = variable)) +
  xlab("Pseudotime") + ylab("") +
  scale_fill_manual(name = "", labels = c("Cell 1", "Cell 2", "Cell 3", "Cell 4"),
                    values = brewer.pal(4, "Set1")) +
  guides(fill = guide_legend(override.aes = list(colour = NULL))) +
  theme(legend.text = element_text(size = 16)) +
  theme(legend.key.height=unit(2, "line"), legend.key.width = unit(2, "line")) +
  cowplot::theme_cowplot()

ggsave(plt, filename = "figs/postuncert/pu_density.png", 
                                 width = 8, height = 2, scale = 1.5)

makeBoxplot <- function(pst) {
  tmcmc <- mcmc(pst)
  hpd75 <- HPDinterval(tmcmc, 0.75)
  hpd95 <- HPDinterval(tmcmc, 0.95)
  
  p <- cbind(hpd75, hpd95)
  p <- data.frame(p)
  names(p) <- c("lower75", "upper75", "lower95", "upper95")
  p$Median <- colMedians(pst)
  p$Cell <- as.factor(rank(p$Median))
  
  
  plt <- ggplot(p) + geom_boxplot(aes(x = Cell, middle = Median, lower = lower75, upper = upper75,
                               ymin = lower95, ymax = upper95), stat = "identity", fill = "darkred", alpha = 0.5) +
      cowplot::theme_cowplot() +
      theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank()) +
       ylab("Pseudotime") # + xlab("Cell") 
  return(plt)
}

post_tracefiles <- paste0("data/", c("trapnell", "burns", "shin"), "_pseudotime_traces.h5")

plts <- lapply(post_tracefiles, function(ptf) {
  pst <- h5read(ptf, "pst")
  makeBoxplot(pst)
})

output_file <- "figs/postuncert/3bcd_post_uncert.png"
pg <- cowplot::plot_grid(plotlist = plts, nrow = 1, labels = c("B","C","D"), rel_widths = c(4, 3, 3))
cowplot::ggsave(pg, filename = output_file, width = 8, height = 2, scale = 1.5)


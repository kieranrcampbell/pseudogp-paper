
library(ggplot2)
library(reshape2)
library(dplyr)
library(RColorBrewer)
library(grid)
library(coda)
setwd("/net/isi-scratch/kieran/GP/pseudogp2/stan/gbio//postuncert")
source("../../gputils//gputils.R")

post_tracefile <- "/net/isi-scratch/kieran/GP/pseudogp2/data/stan_traces_for_gbio.h5"

pst <- h5read(post_tracefile, "pst")

pmean <- colMeans(pst)

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
  theme(legend.key.height=unit(2, "line"), legend.key.width = unit(2, "line"))

ggsave(plt, filename = "pu_density.png", width = 8, height = 2, scale = 1.5)


tmcmc <- mcmc(pst[,ind])
hpd <- HPDinterval(tmcmc, 0.95)

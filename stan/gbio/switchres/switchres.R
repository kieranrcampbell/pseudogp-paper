
# Analysis of the resolution of which a gene can be said to
# switch on or switch off given posterior uncertainty in pseudotime

library(rhdf5)
library(sctools)
library(coda)
library(MCMCglmm)
library(dplyr)
library(GGally)

base_dir <- "~/mount/"
setwd(paste0(base_dir, "GP/pseudogp2/stan/gbio/switchres"))
source("../../diffexpr/monocle/prep_data.R")

sce <- load_data(base_dir)

pst_trace_file <- paste0(base_dir, "GP/pseudogp2/data/stan_traces_for_gbio.h5")
pst <- h5read(pst_trace_file, "pst")
tmap <- posterior.mode(mcmc(pst))

sce$pseudotime <- tmap


# First do differential expression test at map ----------------------------
diff_expr <- testDE(sce)
diff_expr <- data.frame(t(diff_expr))
diff_expr$gene <- fData(sce)$gene_short_name

## now subset down to things that are differentially expressed and have a t0 between 0 and 1
diff_expr$qval <- p.adjust(diff_expr$pval, method = "BH")
de_filter <- filter(diff_expr, qval < 0.05, t0 > 0.05, t0 < 0.95)

ggpairs(de_filter)  #check for spurious correlations
qplot(de_filter$t0, geom = 'density')

## there's an interesting clustering of genes that show strong switch like behaviour and not:
abslog <- function(x) sign(x) * log10(abs(x))
de_filter$kl <- abslog(de_filter$k)
de_filter <- mutate(de_filter, is_switch = abs(kl) > 2.5)
ggplot(de_filter) + geom_point(aes(x = t0, y = kl, color = is_switch)) + xlab("Activation time") + ylab("Activation strength (log10)")

de_strong <- filter(de_filter, abs(kl) > 2.5)
strong_inds <- match(de_strong$gene, fData(sce)$gene_short_name)
sce_strong <- sce[strong_inds,]

strong_fits <- fitModel(sce_strong)

plotlist <- list()
for(i in 1:12) {
  plotlist[[i]] <- norm_plot_model(list(par = strong_fits[,i]), exprs(sce_strong)[i,], sce_strong$pseudotime)
}
plotlist <- plotlist[!sapply(plotlist, is.null)]

cowplot::plot_grid(plotlist = plotlist[!is.null(plotlist)])

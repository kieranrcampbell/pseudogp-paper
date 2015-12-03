
# Analysis of the resolution of which a gene can be said to
# switch on or switch off given posterior uncertainty in pseudotime

library(rhdf5)
library(switchde)
library(coda)
library(MCMCglmm)
library(dplyr)
library(GGally)
library(matrixStats)
library(grid)

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
de_filter  <- de_filter %>% group_by(gene) %>% filter(row_number() == 1) # strip out a duplicate

ggpairs(de_filter)  #check for spurious correlations
qplot(de_filter$t0, geom = 'density')

absqrt <- function(x) sign(x) * sqrt(abs(x))
de_filter$ks <- absqrt(de_filter$k)
de_filter<- mutate(de_filter, sswitch = abs(ks) > 10)

filter_inds <- match(de_filter$gene, fData(sce)$gene_short_name)
sce_filter <- sce[unique(filter_inds),]


to_sample <- sample(nrow(pst), 100)
pstpars <- apply(pst[to_sample,], 1, function(t) {
  sce_strong$pseudotime <- t
  fitModel(sce_strong)[3,]
})

pmodes <- posterior.mode(mcmc(t(pstpars)))
pmad <- apply(pstpars, 1, mad)
qplot(pmad, geom='density')



# switch genes with error bars for *all* filtered genes --------------------------------------------
set.seed(123)
ebar_sample <- sample(nrow(pst), 100)

k_samples <- apply(pst[ebar_sample,], 1, function(t) {
  sce_filter$pseudotime <- t
  switchde::fitModel(sce_filter)[2,]
})

intervals <- HPDinterval(mcmc(t(k_samples)))
intervals <- absqrt(intervals)
idf <- data.frame(intervals, t0 = de_filter$t0)

med_act <- rowMedians(k_samples)
de_filter <- ungroup(de_filter)
de_filter$med_act <- absqrt(med_act)

de_filter <- mutate(de_filter, sswitch = abs(med_act) > 10)
actplt <- ggplot() +
  geom_point(data = de_filter, aes(x = t0, y = med_act), size = 4, alpha = 0.5, shape=21, color='white', fill='black') +
  xlab("Activation time") + ylab(expression(sqrt("Median activation strength"))) + theme_bw() +
  geom_errorbar(data = idf, aes(x = t0, ymin = lower, ymax = upper), alpha = 0.2)
ggsave("5_switchres_a.png", actplt, width=6, height=4)

bottom_switch <- order(abs(de_filter$med_act))[1:5]
top_switch <- order(abs(de_filter$med_act))[nrow(de_filter):(nrow(de_filter)-4)]
de_small <- de_filter[c(top_switch, bottom_switch),]

inds <- match(de_small$gene, fData(sce_filter)$gene_short_name)
scecc <- sce_filter[inds,]

de_filter <- mutate(de_filter, abs_med_act = (abs(med_act)))
de_filter <- mutate(de_filter, mean_exprs = fData(sce_filter)$mean_exprs)

ggpairs(select(de_filter, L, abs_med_act, mean_exprs))

plots <- list()
for(j in 1:nrow(de_small)) {
  current_gene <- inds[j]

  outlier_samples <- apply(pst[ebar_sample,], 1, function(t) {
    sce_filter$pseudotime <- t
    # print(current_gene)
    list(pars = as.numeric(switchde::fitModel(sce_filter[current_gene,])), t = t)
  })

  x <- exprs(sce_filter)[current_gene,]

  plt <- ggplot()
  for(i in 1:length(outlier_samples)) {
    os <- outlier_samples[[i]]
    pars <- os$pars ; t <- os$t
    d <- data.frame(mu = switchde::calc_mu(pars, t), t = t)
    plt <- plt + geom_line(data = d, aes(x = t, y = mu), color = 'red', alpha = 0.3, size=1.2) +
      theme_bw() +
      theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
            plot.margin = unit(c(0.8, 0.1, 0.2, 0.1), "cm"))
  }
  plt <- plt + geom_point(data = data.frame(x, tmap), aes(x = tmap, y = x),
                            size = 3, shape=21, fill='black', color='white', alpha=0.8)
  plots[[j]] <- plt
}
cplt <- cowplot::plot_grid(plotlist=plots, nrow=2, labels=de_small$gene)
ggsave("5_switchres_b.png", cplt, width=6, height=2.5, scale = 2)



# genes that turn on at same time -----------------------------------------

dt <- dist(de_filter$t0)
image(as.matrix(dt))
library(gplots)
heatmap.2(1 - as.matrix(dt), dendrogram="none", trace="none")


# quick volcano plot ------------------------------------------------------

dev <- mutate(de_filter, logqval = -log10(qval))
volcano_plt <- ggplot(dev) + geom_point(aes(x = med_act, y = logqval), alpha = 0.8) +
  theme_bw() + xlab("Median activation strength") + ylab("-log10 Q value")
ggsave("../../supplementary/switch_volcano.png", volcano_plt, width=6, height=4.5)

  # what’s the pseudotemporal switching resolution? -------------------------

t_samples <- apply(pst[ebar_sample,], 1, function(t) {
  sce_filter$pseudotime <- t
  switchde::fitModel(sce_filter)[3,]
})

tmad <- rowMads(t_samples)
tmadmedian <- median(tmad)

madplt <- ggplot(data.frame(tmad)) + geom_density(aes(x = tmad), fill = 'darkred', alpha=0.5) +
  theme_bw() +
  geom_vline(xintercept = tmadmedian, linetype = 2, size = 1.1) + xlim(0, 1) +
  xlab("Median absolute deviation of activation time") + ylab("Density")

ggsave("5c_switchres.png", madplt, width=3, height=2.5, scale = 1.6)

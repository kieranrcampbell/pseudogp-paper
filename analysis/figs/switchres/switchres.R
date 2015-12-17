
# Analysis of the resolution of which a gene can be said to
# switch on or switch off given posterior uncertainty in pseudotime

library(rhdf5)
library(switchde)
library(coda)
library(MCMCglmm)
library(dplyr)
library(matrixStats)
library(grid)
library(scater)


makeSwitchPlots <- function(sce, pst) {
  tmap <- posterior.mode(mcmc(pst))
  
  sce$pseudotime <- tmap
  
  
  # First do differential expression test at map ----------------------------
  diff_expr <- testDE(sce)
  diff_expr <- data.frame(t(diff_expr))
  diff_expr$gene <- featureNames(sce)
  if("gene_short_name" %in% names(fData(sce))) diff_expr$gene <- fData(sce)$gene_short_name
  fData(sce)$gene_short_name <- diff_expr$gene
  
  ## now subset down to things that are differentially expressed and have a t0 between 0 and 1
  diff_expr$qval <- p.adjust(diff_expr$pval, method = "BH")
  de_filter <- filter(diff_expr, qval < 0.05, t0 > 0.05, t0 < 0.95)
  de_filter  <- de_filter %>% group_by(gene) %>% filter(row_number() == 1) # strip out a duplicate
  
  # ggpairs(de_filter)  #check for spurious correlations
  # qplot(de_filter$t0, geom = 'density')
  
  absqrt <- function(x) sign(x) * sqrt(abs(x))
  de_filter$ks <- absqrt(de_filter$k)
  de_filter<- mutate(de_filter, sswitch = abs(ks) > 10)
  
  filter_inds <- match(de_filter$gene, fData(sce)$gene_short_name)
  sce_filter <- sce[unique(filter_inds),]

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
  #ggsave("5_switchres_a.png", actplt, width=6, height=4)
  
  bottom_switch <- order(abs(de_filter$med_act))[1:5]
  top_switch <- order(abs(de_filter$med_act))[nrow(de_filter):(nrow(de_filter)-4)]
  de_small <- de_filter[c(top_switch, bottom_switch),]
  
  inds <- match(de_small$gene, fData(sce_filter)$gene_short_name)
  scecc <- sce_filter[inds,]
  
  de_filter <- mutate(de_filter, abs_med_act = (abs(med_act)))
  de_filter <- mutate(de_filter, mean_exprs = fData(sce_filter)$mean_exprs)
  
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

  return(list(actplt, cplt, madplt))
}



base_dir <- "/net/isi-scratch/kieran/"

load(file.path(base_dir, "pseudogp-paper/data/sce_monocle.Rdata"))
sce_monocle <- sce_23

load(file.path(base_dir, "pseudogp-paper/data/sce_ear.Rdata"))
sce_ear <- sct

load(file.path(base_dir, "pseudogp-paper/data/sce_waterfall.Rdata"))
sce_waterfall <- sce

sces <- list(monocle = sce_monocle, ear = sce_ear, waterfall = sce_waterfall)

post_tracefiles <- paste0(base_dir, "pseudogp-paper/data/",
                          c("monocle_stan_traces.h5", "ear_stan_traces.h5", "waterfall_stan_traces.h5"))

psts <- lapply(post_tracefiles, function(ptf) {
  h5read(ptf, "pst")
})

all_plts <- lapply(1:3, function(i) {
  sce <- sces[[i]]
  n_cells_exprs <- rowSums(exprs(sce) > sce@lowerDetectionLimit)
  genes_to_use <- n_cells_exprs > (0.1 * ncol(sce)) # select genes expressed in at least 10% of cells
  sce <- sce[genes_to_use,]
  makeSwitchPlots(sce, psts[[i]])
})

ns <- c("monocle","ear","waterfall")
setwd(file.path(base_dir, "pseudogp-paper/analysis/figs/switchres/"))
for(i in 1:3) {
  base_name <- paste0(ns[i], "_5_switchres")
  plts <- all_plts[[i]]
  ggsave(paste0(base_name, "_a.png"), actplt, width=6, height=4)
  ggsave(paste0(base_name, "_b.png"), cplt, width=6, height=2.5, scale = 2)
  ggsave(paste0(base_name, "_c.png"), madplt, width=3, height=2.5, scale = 1.6)
}


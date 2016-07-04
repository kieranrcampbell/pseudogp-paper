
# Differential expression comparison for uncertainty in pseudotime
# kieran.campbell@sjc.ox.ac.uk

library(ggplot2)
library(readr)
library(magrittr)
library(magrittr)
library(ggplot2)
library(dplyr)
library(tidyr)

alpha <- 0.05 # significance level
# study <- "shin"

# args <- commandArgs(trailingOnly = TRUE)
# study <- args[1]
# csv_file <- args[2]
# stopifnot(study %in% c("trapnell", "burns", "shin"))


read_studies <- function(study) {
  agg_pvals <- paste0("data/diffexpr/agg_pvals/", study, ".csv")
  map_pvals <- paste0("data/diffexpr/map/", study, ".csv")
  
  ap <- read_csv(agg_pvals)
  mapp <- read_csv(map_pvals)
  
  prop_sig <- ap %>% group_by(gene) %>%
    summarise(prop_sig = mean(q_val < alpha))
  
  sig_df <- inner_join(mapp, prop_sig, by = "gene")
  
  designate <- function(qval, prop_sig) {
    if(qval < alpha & prop_sig >= (1 - alpha)) {
      return("True positive")
    } else if(qval < alpha & prop_sig < (1 - alpha)) {
      return("False positive")
    } else if(qval > alpha & prop_sig < (1 - alpha)) {
      return("True negative") 
    } else {
      return("False negative")
    }
  }
  
  discovery_types <- sig_df %>%
    select(q_val, prop_sig) %>%
    apply(1, function(row) designate(row[1], row[2]))
  
  sig_df %<>% mutate(discovery_type = discovery_types)
  sig_df %<>% mutate(percent_sig = 100 * prop_sig)
  sig_df %<>% mutate(study = study)
  return(sig_df)
}

studies <- c("trapnell", "shin", "burns")
study_list <- lapply(studies, read_studies)
sig_df <- do.call(rbind, study_list)
write_csv(sig_df, "data/diffexpr/all_pvals.csv")

sig_df <- read_csv("data/diffexpr/all_pvals.csv")
fd_df <- sig_df %>%
  group_by(study) %>%
  summarise(n_fp = sum(discovery_type == "False positive"),
            n_tp = sum(discovery_type == "True positive"),
            n_tn = sum(discovery_type == "True negative")) %>%
  mutate(FDR = n_fp / (n_fp + n_tn))


# allgene_plt <- ggplot(sig_df, aes(x = percent_sig, y = q_val, fill = discovery_type)) +
#   geom_point(alpha = 0.5, shape = 21, color = 'black', size = 2.3) +
#   cowplot::theme_cowplot() +
#   scale_fill_brewer(name = element_blank(), palette = "Set1") +
#   xlab("% MCMC traces significant") +
#   ylab("FDR-adjusted MAP p-value") +
#   theme(legend.position = c(0.7,0.8)) +
#   geom_hline(yintercept = alpha, colour = "grey", linetype = 2) +
#   geom_vline(xintercept = 100 * (1 - alpha), colour = "grey", linetype = 2) 

fd_df$study <- plyr::mapvalues(fd_df$study, from = c("burns", "shin", "trapnell"),
                           to = c("Burns", "Shin", "Trapnell"))


f <- fd_df %>% 
  gather(measure, number, -FDR, -study) %>%
  filter(measure != "n_tn") %>% 
  arrange(desc(measure))

f$measure <- plyr::mapvalues(f$measure, from = c("n_tp", "n_fp"),
                             to = c("True positive", "False positive"))
# f$measure <- factor(f$measure, levels = c("True postive", "False positive"))


de_gene_nums <- ggplot(f, aes(x = study, y = number, fill = measure)) + 
  geom_bar(stat = "identity") +
  cowplot::theme_cowplot() +
  scale_fill_brewer(name = element_blank(), palette = "Set1") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.8, size = 10, hjust = 0.8)) +
  ylab("Number of DE genes") 

# fd_df <- mutate(fd_df, fdr_pct = paste0(round(100 * FDR, digits = 1), "%"))

# fdr_plt <- ggplot(fd_df, aes(x = study, y = FDR, fill = study)) + geom_bar(stat = "identity") +
#   cowplot::theme_cowplot() +
#   xlab("Study") + ylab("False discovery rate") +
#   scale_fill_brewer(name = element_blank(), palette = "Set1") +
#   theme(legend.position = "none") +
#   geom_text(aes(label = fdr_pct), vjust = -0.25)

# ggsave(plot = allgene_plt, filename = "figs/diffexpr/all_genes.png", width = 5.5, height = 5)
ggsave(plot = de_gene_nums, filename = "figs/diffexpr/de_gene_nums.png", width = 6, height = 4)
# ggsave(plot = fdr_plt, filename = "figs/diffexpr/fdr_barplot.png", width = 5.5, height = 5)


# Figures 4 A and B - ITGAE and ID1 examples ------------------------------

# Helper functions


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




## Analysis of p-values

library(readr)
library(ggplot2)
library(cowplot)
library(modeest)
library(matrixStats)
library(reshape2)
library(scater)
library(embeddr)
library(dplyr)

source("prep_data.R")
source("/net/isi-scratch/kieran/switch/sctools/R/normal_switch.R")

d <- "/net/isi-scratch/kieran/GP/pseudogp2/data/ldl01/"
files <- dir(d)
files <- files[grep("pval", files)]
ss_files <- files[grep("ss", files)]
switch_files <- files[grep("switch", files)]

ss_pvals <- switch_pvals <- list()

for(f in ss_files) ss_pvals <- c(ss_pvals, read_csv(paste0(d, f)))
for(f in switch_files) switch_pvals <- c(switch_pvals, read_csv(paste0(d, f)))

ss_p <- do.call("cbind", ss_pvals)
switch_p <- do.call("cbind", switch_pvals)



## normal_switch gives a p-value of -1 if the alternative model couldn't be fit
## and -1 if the null model couldn't be fit
optfail_per_gene <- apply(switch_p, 1, function(x) sum(x == -1))
nullfail_per_gene <- apply(switch_p, 1, function(x) sum(x == -2)) ## this should be none

print(sum(optfail_per_gene))
print(sum(nullfail_per_gene))

switch_p_copy <- switch_p
switch_p[switch_p == -1] <- 1

ss_q <- apply(ss_p, 2, p.adjust, "BH")
switch_q <- apply(switch_p, 2, p.adjust, "BH")

ss_is_sig <- ss_q < 0.05
switch_is_sig <- switch_q < 0.05

ss_med_q <- rowMedians(ss_q)
sw_med_q <- rowMedians(switch_q)

ss_prop_sig <- rowSums(ss_is_sig) / ncol(ss_p)
switch_prop_sig <- rowSums(switch_is_sig) / ncol(switch_p)

df <- data.frame(ss = ss_prop_sig, sw = switch_prop_sig)
scatter_plt <- ggplot(df) + 
  geom_point(aes(x = ss, y = sw), alpha = 0.3, size = 1) +
  xlab("Smoothing splines") + ylab("Switch-like")

names(df) <- c("Smoothing splines", "Switch-like")
dfm <- melt(df, variable.name = "test", value.name = "prop_sig")
hist_plt <- ggplot(dfm) + geom_histogram(aes(x = prop_sig)) + facet_wrap(~ test)  + xlab("Proportion significant")


ss_grey <- sum(ss_prop_sig > 0 & ss_prop_sig < 1)
sw_grey <- sum(switch_prop_sig > 0 & switch_prop_sig < 1)

df_grey <- data.frame(test = factor(c("Smoothing splines", "Switch")), count = c(ss_grey, sw_grey))

bar_plt <- ggplot(df_grey) + geom_bar(aes(x = test, y = count), stat="identity") +
  xlab("Differential expression test") + ylab("")

grid_plt <- plot_grid(hist_plt, scatter_plt, ncol = 1, labels = c("A", "B", "C", "D"))
cowplot::ggsave("/net/isi-scratch/kieran/GP/pseudogp2/diffexpr/pvalsplot.png", 
                grid_plt, scale = 1.5, width = 4, height = 6)


# Compare representative samples of splines vs switch ------
names(df) <- c('ss', 'sw')
gene_type <- rep('agree', nrow(df))
gene_type[df$sw > 0.8 & df$ss < 0.2] <- "switch"
gene_type[df$ss > 0.99 & df$sw == 0] <- "smoothing"
table(gene_type)
df$gene_type <- gene_type
ggplot(df) + 
  geom_point(aes(x = ss, y = sw, colour = gene_type), alpha = 0.3, size = 3) +
  xlab("Smoothing splines") + ylab("Switch-like")

get_map <- function(x) mlv(x,method = "HSM")$M

data <- load_data()
sce <- data$sce ; pst <- data$pst

pst_map <- apply(pst, 2, get_map)
sce$pseudotime <- pst_map

## let's look at those funky genes - first good switch, bad smoothing
good_switch <- gene_type == "switch"
# ss_test <- pseudotime_test(sce[good_switch,], n_cores = 1)
ss_models <- plot_pseudotime_model(sce[good_switch,], n_cores = 1, ncol = 3, 
                                   facet_wrap_scale = "free", mask_min_expr = TRUE)

switch_plots <- lapply(which(good_switch), function(j) {
  x <- as.vector(exprs(sce)[j,])
  model <- fit_alt_model(x, pseudotime(sce))
  model$type <- "sigmoid"
  return( plot_model(model, x, pseudotime(sce)) )
})

plot_grid(ss_models, plot_grid(plotlist = switch_plots, ncol = 3), ncol = 2)

## now ss model, bad switching
good_ss <- gene_type == "smoothing"
set.seed(123)
good_ss <- sample(which(good_ss), sum(good_switch)) # pick a number of equal size to good switch

ss_models2 <- plot_pseudotime_model(sce[good_ss,], n_cores = 1, ncol = 3, 
                                    facet_wrap_scale = "free", mask_min_expr = TRUE)

switch_plots2 <- lapply(good_ss, function(j) {
  x <- as.vector(exprs(sce)[j,])
  model <- fit_alt_model(x, pseudotime(sce))
  model$type <- "sigmoid"
  return( plot_model(model, x, pseudotime(sce)) )
})

plot_grid(ss_models2, plot_grid(plotlist = switch_plots2, ncol = 3), ncol = 2)


# Time for analysis of FDR -----------

de_test <- pseudotime_test(sce, n_cores = 1)
# compare to ss_prop_sig
dfc <- data.frame(psig = ss_prop_sig, qval = de_test$q_val, med_q_val = ss_med_q)
dfc <- mutate(dfc, is_sig = qval < 0.05)

plt_ss_prop <- ggplot(dfc) + geom_point(aes(x = psig, y = qval, colour = is_sig), size = 2, alpha = 0.5) +
  ggthemes::scale_colour_economist() + 
  geom_hline(yintercept = 0.05, linetype = 1, colour = 'darkred') + 
  ylab("Q-val for MAP pseudotime estimate") + xlab(" ") +
  theme(legend.position = "none") + 
  geom_rug(data = filter(dfc, is_sig), aes(x = psig, y = qval), sides = "b", alpha = 0.1, colour = "darkred")

plt_ss_median <- ggplot(dfc) + geom_point(aes(x = med_q_val, y = qval, colour = is_sig), size = 2, alpha = 0.5) +
  ggthemes::scale_colour_economist() + 
  geom_hline(yintercept = 0.05, linetype = 1, colour = 'darkred') + 
  geom_vline(xintercept = 0.05, linetype = 2, colour = 'darkred') +
  ylab(" ") + xlab(" ") +
  theme(legend.position = "none") + 
  geom_rug(data = filter(dfc, is_sig), aes(x = med_q_val, y = qval), sides = "b", alpha = 0.1, colour = "darkred")


cowplot::ggsave(filename = "fdr.pdf")

## plot for Chris with density
cplt <- ggplot(dfc) + geom_point(aes(x = psig, y = qval, colour = is_sig), size = 2, alpha = 0.5) +
  ggthemes::scale_colour_economist() + 
  geom_hline(yintercept = 0.05, linetype = 1, colour = 'darkred') + 
  ylab("Q-val for MAP pseudotime estimate") + 
  ggplot2::theme_bw() +
  theme(legend.position = "none",          
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        #axis.text.x = element_blank(),
        axis.text.y = element_blank(), 
        plot.margin = grid::unit(c(3,-5.5,4,3), "mm"))

hist_top <- ggplot(filter(dfc, is_sig)) + geom_density(aes(x = psig), fill = "darkgrey") + 
  ggplot2::theme_bw() +   
  theme(legend.position = "none",          
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        #axis.text.x = element_blank(),
        axis.text.y = element_blank(), 
        plot.margin = grid::unit(c(3,-5.5,4,3), "mm")) +
  scale_x_continuous(limits = c(0, 1))
pdf("cplt.pdf")
gridExtra::grid.arrange(cplt, hist_top, ncol = 1, heights = c(4,1))
dev.off()

## repeat the same for switch-like
switch_pvals <- apply(exprs(sce), 1, function(x) {
  diff_expr_test(x, pseudotime(sce))[1]
})
stopifnot(all(switch_pvals >= 0))
switch_qvals <- p.adjust(switch_pvals, method = "BH")
dfs <- data.frame(psig = switch_prop_sig, qval = switch_qvals, med_q_val = sw_med_q)
dfs <- mutate(dfs, is_sig = qval < 0.05)

plt_sw_prop <- ggplot(dfs) + geom_point(aes(x = psig, y = qval, colour = is_sig), size = 2, alpha = 0.5) +
  ggthemes::scale_colour_economist() + 
  geom_hline(yintercept = 0.05, linetype = 1, colour = 'darkred') + 
  ylab("Q-val for MAP pseudotime estimate") + xlab("Proportion significant (FDR 5%)") +
  theme(legend.position = "none") + 
  geom_rug(data = filter(dfs, is_sig), aes(x = psig, y = qval), sides = "b", alpha = 0.1, colour = "darkred")

plt_sw_median <- ggplot(dfs) + geom_point(aes(x = med_q_val, y = qval, colour = is_sig), size = 2, alpha = 0.5) +
  ggthemes::scale_colour_economist() + 
  geom_hline(yintercept = 0.05, linetype = 1, colour = 'darkred') + 
  geom_vline(xintercept = 0.05, linetype = 2, colour = 'darkred') +
  ylab(" ") + xlab("Median Q-value") +
  theme(legend.position = "none") + 
  geom_rug(data = filter(dfs, is_sig), aes(x = med_q_val, y = qval), sides = "b", alpha = 0.1, colour = "darkred")

all_plot <- plot_grid(plt_ss_prop, plt_ss_median, plt_sw_prop, plt_sw_median, ncol = 2)
cowplot::ggsave("all.pdf", plot = all_plot)

## boxplot time
ss_bxplt <- select(filter(dfc, is_sig), psig, med_q_val)
sw_bxplt <- select(filter(dfs, is_sig), psig, med_q_val)

#sw_bxplt <- rbind(sw_bxplt, matrix(NA, nrow = nrow(ss_bxplt) - nrow(sw_bxplt), ncol = 2))
d <- rbind(ss_bxplt, sw_bxplt)
d$type <- as.factor(c(rep("Smoothing spline", nrow(ss_bxplt)), rep("Switch-like", nrow(sw_bxplt))))
dm <- melt(d, id.vars = "type", variable.name = "metric")
dm$metric <- plyr::mapvalues(dm$metric, from = c("psig", "med_q_val"), to = c("Proportion significant", "Median Q value"))

# df_bxplt <- data.frame("Smoothing splines" = ss_bxplt, "Switch-like" = sw_bxplt)
# dfbm <- melt(df_bxplt, variable.name = "test", value.name = "prop_sig")

ggplot(dm, aes(y = value, x = type)) + geom_boxplot() + 
  geom_point(position = position_jitter(width = 0.22), alpha = 0.1) +
  xlab("Differential gene test") + ylab("") +
  facet_wrap(~ metric) + cowplot::theme_cowplot() 


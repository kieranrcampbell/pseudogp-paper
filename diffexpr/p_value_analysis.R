
## Analysis of p-values

library(readr)
library(ggplot2)
library(cowplot)
library(modeest)

source("prep_data.R")

d <- "/net/isi-scratch/kieran/GP/gpseudotime/data/decsv/"
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

switch_p_copy <- switch_p
switch_p[switch_p == -1] <- 1

ss_is_sig <- ss_p < 0.05
switch_is_sig <- switch_p < 0.05

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
names(df) <- c('sw', 'ss')
gene_type <- rep('agree', nrow(df))
gene_type[df$sw > 0.9 & df$ss < 0.1] <- "switch"
gene_type[df$ss > 0.995 & df$sw == 0.0] <- "smoothing"
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
ss_models <- plot_pseudotime_model(sce[good_switch,], n_cores = 1, ncol = 3, facet_wrap_scale = "free")

switch_plots <- lapply(which(good_switch), function(j) {
  x <- as.vector(exprs(sce)[j,])
  model <- fit_alt_model(x, pseudotime(sce))
  model$type <- "sigmoid"
  return( plot_model(model, x, pseudotime(sce)) )
})

plot_grid(ss_models, plot_grid(plotlist = switch_plots, ncol = 3), ncol = 2)

## now ss model, bad switching
good_ss <- gene_type == "smoothing"
good_ss <- sample(which(good_ss), sum(good_switch)) # pick a number of equal size to good switch

ss_models2 <- plot_pseudotime_model(sce[good_ss,], n_cores = 1, ncol = 3, facet_wrap_scale = "free")

switch_plots2 <- lapply(good_ss, function(j) {
  x <- as.vector(exprs(sce)[j,])
  model <- fit_alt_model(x, pseudotime(sce))
  model$type <- "sigmoid"
  return( plot_model(model, x, pseudotime(sce)) )
})

plot_grid(ss_models2, plot_grid(plotlist = switch_plots2, ncol = 3), ncol = 2)





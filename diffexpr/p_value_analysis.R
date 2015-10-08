
## Analysis of p-values

library(readr)
library(ggplot2)
library(cowplot)


d <- "~/mount/GP/gpseudotime/data/decsv/"
files <- dir(d)
files <- files[grep("pval", files)]
ss_files <- files[grep("ss", files)]
switch_files <- files[grep("switch", files)]

ss_pvals <- switch_pvals <- list()

for(f in ss_files) ss_pvals <- c(ss_pvals, read_csv(paste0(d, f)))
for(f in switch_files) switch_pvals <- c(switch_pvals, read_csv(paste0(d, f)))

ss_p <- do.call("cbind", ss_pvals)
switch_p <- do.call("cbind", switch_pvals)

ss_is_sig <- ss_p < 0.05
switch_is_sig <- switch_p < 0.05

ss_prop_sig <- rowSums(ss_is_sig) / ncol(ss_p)
switch_prop_sig <- rowSums(switch_is_sig) / ncol(switch_p)

df <- data.frame(ss = ss_prop_sig, sw = switch_prop_sig)

ss_plt <- ggplot(df) + geom_histogram(aes(x = ss)) + xlab("Proportion significant")
sw_plt <- ggplot(df) + geom_histogram(aes(x = sw)) + xlab("Proportion significant")
scatter_plt <- ggplot(df) + geom_point(aes(x = ss, y = sw), alpha = 0.3, size = 1) +
  xlab("Smoothing splines") + ylab("Switch-like")

ss_grey <- sum(ss_prop_sig > 0 & ss_prop_sig < 1)
sw_grey <- sum(switch_prop_sig > 0 & switch_prop_sig < 1)

df_grey <- data.frame(test = factor(c("Smoothing splines", "Switch")), count = c(ss_grey, sw_grey))

bar_plt <- ggplot(df_grey) + geom_bar(aes(x = test, y = count), stat="identity") +
  xlab("Differential expression test") + ylab("")

grid_plt <- plot_grid(ss_plt, sw_plt, scatter_plt, bar_plt, ncol = 2, labels = c("A", "B", "C", "D"))
cowplot::ggsave("~/oxford/GP/pseudogp2/diffexpr/pvalsplot.png", grid_plt, scale = 2)


# Bit more in the grey zone -----------------------------------------------

grey_genes <- ss_prop_sig > 0 & ss_prop_sig < 1
grey_p <- ss_p[grey_genes,]


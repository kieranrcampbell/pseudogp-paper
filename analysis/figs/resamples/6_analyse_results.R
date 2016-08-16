library(dplyr)
library(ggplot2)
library(cowplot)

alpha <- 0.05

rdir <- "data/resamples/trace_diffexpr/"
files <- dir(rdir)

df_list <- lapply(files, function(f) {
  split <- strsplit(f, "_")[[1]]
  pca_rep <- split[2]
  df <- read_csv(file.path(rdir, f))
  df <- dplyr::mutate(df, pca_rep = pca_rep)
  return(df)
})

rdf <- bind_rows(df_list)
rm(df_list)

psig_df <- rdf %>%
  group_by(gene, pca_rep) %>%
  summarise(prop_sig = mean(q_val < alpha))

rm(rdf)
gc()

psig_df %<>% mutate(robust_sig = prop_sig > (1 - alpha))

p2 <- psig_df %>% 
  group_by(gene) %>% 
  summarise(mean_prop_sig = mean(prop_sig),
            prop_pca_rsig = mean(robust_sig == TRUE))



plt_all_traces <- ggplot(psig_df, aes(x = 100 * prop_sig)) + geom_histogram(aes(y = ..count../sum(..count..))) +
  xlab("% DE tests") +
  ylab("Frequency")

print(mean(psig_df$robust_sig))


plt_gene <-   ggplot(p2, aes(x = 100 * prop_pca_rsig)) + 
  geom_histogram(aes(y = cumsum(..count../sum(..count..)))) +
  xlab("% PCA subsamples") +
  ylab("Cumulative frequency")

prop_pca <- p2$prop_pca_rsig
print(summary(prop_pca))
print(mean(prop_pca == 1))
print(mean(prop_pca > 0.8))


by_pca <- psig_df %>%
  group_by(pca_rep) %>%
  summarise(mean_per_pca = mean(robust_sig))

plt_pca <- ggplot(by_pca, aes(x = mean_per_pca)) + geom_histogram() +
  xlab("Proportion of genes robustly DE \n in each PCA subsample") +
  ylab("Count")

sig_genes <- by_pca$mean_per_pca
print(summary(sig_genes))


plot_grid(plt_all_traces, plt_gene, nrow = 1, labels = "AUTO")

ggsave("figs/resample_results.png", width = 7, height = 3.5)

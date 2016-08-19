library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(MCMCglmm)
library(coda)
library(cowplot)
library(goseq)
library(magrittr)

alpha <- 0.05

bootstrap_dir <- "data/bootstrap/bootstrapped_de/"
gplvm_dir <- "data/bootstrap/gplvm_de"

bootstrap_files <- dir(bootstrap_dir)
bootstrap_files <- bootstrap_files[grep("de", bootstrap_files)]
bootstrap_files <- bootstrap_files[grep("0", bootstrap_files)]


bootstrap_dfs <- lapply(bootstrap_files, function(bs) {
  parsed <- strsplit(bs, "_")[[1]]
  prop_sampled <- as.numeric(parsed[2])
  boostrap_no <- gsub(".csv", "", parsed[3], fixed = TRUE)
  dx <- read_csv(file.path(bootstrap_dir, bs))
  dx <- mutate(dx, Bootstrap = bootstrap_no, prop_sampled)
  dx
})

bootstrap_df <- bind_rows(bootstrap_dfs)
rm(bootstrap_dfs)

gplvm_files <- dir(gplvm_dir)
gplvm_files <- gplvm_files[grep("de", gplvm_files)]

gplvm_dfs <- lapply(gplvm_files, function(bs) {
  gplvm_no <- as.integer(gsub("de_", "", gsub(".csv", "", bs, fixed = TRUE)))
  dx <- read_csv(file.path(gplvm_dir, bs))
  dx <- mutate(dx, gplvm = gplvm_no)
  dx
})

gplvm_df <- bind_rows(gplvm_dfs)
rm(gplvm_dfs)

bootstrap_psig <- group_by(bootstrap_df, gene, prop_sampled) %>% 
  summarise(psig = mean(q_val < alpha), mean_qval = mean(q_val))

gplvm_psig <- group_by(gplvm_df, gene) %>% 
  summarise(psig = mean(q_val < alpha), mean_qval = mean(q_val))

df <- inner_join(bootstrap_psig, gplvm_psig, by = "gene", suffix = c("_bootstrap", "_gplvm"))

df$prop_sampled <- plyr::mapvalues(df$prop_sampled, from = c(0.5, 0.7, 0.9),
                                   to = c("50%", "70%", "90%"))

pltA <- ggplot(df, aes(x = psig_bootstrap, y = psig_gplvm)) + geom_point(alpha = 0.4) +
  xlab("Proportion significant bootstrap") +
  ylab("Proportion significant \nGPLVM") +
  facet_wrap(~ prop_sampled)

ggsave("~/Desktop/reviewer1.png", width = 10, height = 4)

pltB <- ggplot(df, aes(x = mean_qval_bootstrap, y = mean_qval_gplvm)) + geom_point(alpha = 0.4) +
  xlab("Mean Q-val bootstrap") +
  ylab("Mean Q-val \nGPLVM") + facet_wrap(~ prop_sampled) +
  stat_function(fun = function(x) x, color = 'red')

ggsave("~/Desktop/reviewer2.png", width = 10, height = 4)


# More plots

disp_df <- filter(df, prop_sampled == 0.9, psig_gplvm < 0.75) %>% 
  mutate(disparity = abs(psig_bootstrap - psig_gplvm)) %>% 
  arrange(desc(disparity))

agree_df <- filter(df, prop_sampled == 0.9, psig_bootstrap > 0.95, psig_gplvm > 0.95) %>% 
  mutate(disparity = abs(psig_bootstrap - psig_gplvm)) %>% 
  arrange(disparity)

bad_ensembl_id <- disp_df$gene[2]
myog_ensembl_id <- featureNames(sce)[grep("MYOG", fData(sce)$gene_short_name)]

filter(agree_df, gene == myog_ensembl_id)


# Time to make some gene plots --------------------------------------------

pseudotime_traces <- read_csv("data/bootstrap/gplvm_pseudotimes.csv")

pst_map <- posterior.mode(mcmc(t(as.matrix(pseudotime_traces))))

gex_bad <- exprs(sce)[bad_ensembl_id, ]
gex_good <- exprs(sce)[myog_ensembl_id, ]

bad_hgnc_symbol <- fData(sce)[bad_ensembl_id, ]$gene_short_name

lower_plt <- data_frame(PSPHP1 = gex_bad, MYOG = gex_good, pst_map) %>% 
  gather(Gene, Expression, -pst_map) %>% 
  ggplot(aes(x = pst_map, y = Expression)) + geom_point() +
  cowplot::theme_cowplot() + stat_smooth(se = FALSE, color = "red") +
  xlab("MAP pseudotime") + ylab("Expression") +
  facet_wrap(~ Gene)


# upper_grid <- plot_grid(pltA, pltB, nrow = 1, labels = "AUTO")
lower_grid <- plot_grid(lower_plt, labels = "C", scale = 0.8)

plot_grid(pltA, pltB, lower_grid, nrow = 3)

ggsave("figs/bootstrap/bootstrap_fig.png", width = 9, height = 8)


# Compare bootstrapped pseudotimes to gplvm pseudotimes -------------------

disp_df_2 <- filter(df, psig_gplvm > 0.95) %>% 
  mutate(disparity = abs(psig_bootstrap - psig_gplvm)) %>% 
  arrange(desc(disparity))

gp_bad_gene <- disp_df_2$gene[6]
gex_gp_bad <- exprs(sce)[gp_bad_gene, ]

data_frame(gex_gp_bad, pst_map) %>% 
  ggplot(aes(x = pst_map, y = gex_gp_bad)) + geom_point() +
  stat_smooth()



# GO time -----------------------------------------------------------------

# UpsetR plot time --------------------------------------------------------

library(UpSetR)

df9 <- filter(df, prop_sampled == 0.9)

ulist = list(bootstrap = filter(df9, psig_bootstrap > 0.95) %>% extract2("gene"),
             gplvm = filter(df9, psig_gplvm > 0.95) %>% extract2("gene"))

png("~/Desktop/upset.png")
upset(fromList(ulist))
dev.off()

library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(MCMCglmm)
library(coda)

alpha <- 0.05

bootstrap_dir <- "data/bootstrap/bootstrapped_de/"
gplvm_dir <- "data/bootstrap/gplvm_de"

bootstrap_files <- dir(bootstrap_dir)
bootstrap_files <- bootstrap_files[grep("de", bootstrap_files)]

bootstrap_dfs <- lapply(bootstrap_files, function(bs) {
  bootstrap_no <- as.integer(gsub("de_", "", gsub(".csv", "", bs, fixed = TRUE)))
  dx <- read_csv(file.path(bootstrap_dir, bs))
  dx <- mutate(dx, Bootstrap = bootstrap_no)
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

bootstrap_psig <- group_by(bootstrap_df, gene) %>% 
  summarise(psig = mean(q_val < alpha), mean_qval = mean(q_val))

gplvm_psig <- group_by(gplvm_df, gene) %>% 
  summarise(psig = mean(q_val < alpha), mean_qval = mean(q_val))

df <- inner_join(bootstrap_psig, gplvm_psig, by = "gene", suffix = c("_bootstrap", "_gplvm"))


ggplot(df, aes(x = psig_bootstrap, y = psig_gplvm)) + geom_point() +
  xlab("Proportion significant using bootstrap") +
  ylab("Proportion significant using GPLVM") + theme_classic()

ggsave("~/Desktop/reviewer1.png")

ggplot(df, aes(x = mean_qval_bootstrap, y = mean_qval_gplvm)) + geom_point() +
  xlab("Mean Q-val bootstrap") +
  ylab("Mean Q-val GPLVM") + theme_classic() +
  stat_function(fun = function(x) x, color = 'red')

ggsave("~/Desktop/reviewer2.png")


# More plots

disp_df <- filter(df, psig_gplvm < 0.75) %>% 
  mutate(disparity = abs(psig_bootstrap - psig_gplvm)) %>% 
  arrange(desc(disparity))

agree_df <- filter(df, psig_bootstrap > 0.95, psig_gplvm > 0.95) %>% 
  mutate(disparity = abs(psig_bootstrap - psig_gplvm)) %>% 
  arrange(disparity)

bad_ensembl_id <- disp_df$gene[1]
myog_ensembl_id <- featureNames(sce)[grep("MYOG", fData(sce)$gene_short_name)]

filter(agree_df, gene == myog_ensembl_id)


# Time to make some gene plots --------------------------------------------



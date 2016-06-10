
# Differential expression comparison for uncertainty in pseudotime
# kieran.campbell@sjc.ox.ac.uk

library(scater)
library(embeddr)
library(readr)
library(MCMCglmm)
library(coda)
library(rhdf5)

alpha <- 0.05

args <- commandArgs(trailingOnly = TRUE)
study <- args[1]
csv_file <- args[2]
stopifnot(study %in% c("trapnell", "burns", "shin"))

agg_pvals <- paste0("data/diffexpr/agg_pvals/", study, ".csv")
map_pvals <- paste0("data/diffexpr/map/", study, ".csv")

ap <- read_csv(agg_pvals)
mapp <- read_csv(map_pvals)

prop_sig <- ap %>% group_by(gene) %>%
  summarise(prop_sig = mean(q_val < alpha))

sig_df <- inner_join(mapp, prop_sig, by = "gene")

ggplot(sig_df, aes(x = q_val, y = prop_sig)) + geom_point()

sig_df %<>% mutate(false_disc = q_val < alpha & prop_sig < (1 - alpha))

print(paste("False discovery rate for study", study, "is", mean(sig_df$false_disc)))

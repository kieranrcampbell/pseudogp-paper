
# Differential expression comparison for uncertainty in pseudotime
# kieran.campbell@sjc.ox.ac.uk

library(ggplot2)
library(readr)
library(magrittr)
library(magrittr)
library(ggplot2)
library(dplyr)

alpha <- 0.05 # significance level

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

num_fp <- sum(sig_df$discovery_type == "False positive")
num_tp <- sum(sig_df$discovery_type == "True positive")

fpr <- num_fp / (num_fp + num_tp)


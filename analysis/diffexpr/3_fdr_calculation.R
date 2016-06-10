
# Differential expression comparison for uncertainty in pseudotime
# kieran.campbell@sjc.ox.ac.uk


library(readr)
library(magrittr)
library(ggplot2)
library(dplyr)

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

print(paste("False discovery ratio for study", study, "is", fpr))

sig_df %<>% mutate(percent_sig = 100 * prop_sig)

ggplot(sig_df, aes(x = percent_sig, y = q_val, fill = discovery_type)) +
  geom_point(alpha = 0.5, shape = 21, color = 'black', size = 2.3) +
  cowplot::theme_cowplot() +
  scale_fill_brewer(name = element_blank(), palette = "Set1") +
  xlab("% MCMC traces significant") +
  ylab("FDR-adjusted MAP p-value") +
  theme(legend.position = c(0.7,0.8)) +
  geom_hline(yintercept = alpha, colour = "grey", linetype = 2) +
  geom_vline(xintercept = 100 * (1 - alpha), colour = "grey", linetype = 2) +
  ggtitle(paste("FDR analysis for", study, "et al."))

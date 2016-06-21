
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


fd_df <- sig_df %>%
  group_by(study) %>%
  summarise(n_fp = sum(discovery_type == "False positive"),
            n_tp = sum(discovery_type == "True positive"),
            n_tn = sum(discovery_type == "True negative")) %>%
  mutate(FDR = n_fp / (n_fp + n_tn))


ggplot(sig_df, aes(x = percent_sig, y = q_val, fill = discovery_type)) +
  geom_point(alpha = 0.5, shape = 21, color = 'black', size = 2.3) +
  cowplot::theme_cowplot() +
  scale_fill_brewer(name = element_blank(), palette = "Set1") +
  xlab("% MCMC traces significant") +
  ylab("FDR-adjusted MAP p-value") +
  theme(legend.position = c(0.7,0.8)) +
  geom_hline(yintercept = alpha, colour = "grey", linetype = 2) +
  geom_vline(xintercept = 100 * (1 - alpha), colour = "grey", linetype = 2) +
  ggtitle(paste("FDR analysis across all studies"))

fd_df$study <- plyr::mapvalues(fd_df$study, from = c("burns", "shin", "trapnell"),
                           to = c("Burns", "Shin", "Trapnell"))


f <- fd_df %>% 
  gather(measure, number, -FDR, -study) %>%
  filter(measure != "n_tn")

f$measure <- plyr::mapvalues(f$measure, from = c("n_fp", "n_tp"),
                             to = c("False positive", "True positive"))

ggplot(f, aes(x = measure, y = number, fill = measure)) + 
  geom_bar(stat = "identity") +
  facet_wrap(~ study) +
  cowplot::theme_cowplot() +
  scale_fill_brewer(name = element_blank(), palette = "Set1") +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.8, size = 10, hjust = 0.8)) +
  ylab("Number")

ggplot(fd_df, aes(x = study, y = FDR)) + geom_bar(stat = "identity") +
  coord_flip() + cowplot::theme_cowplot() +
  xlab("Study") + ylab("False discovery rate")




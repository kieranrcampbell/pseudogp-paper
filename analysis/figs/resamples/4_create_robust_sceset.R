library(scater)
library(readr)
library(dplyr)
library(magrittr)

alpha <- 0.05

de_dir <- "data/resamples/all_cells_diffexpr/"

de_files <- dir(de_dir)

dfs <- lapply(de_files, function(f) {
  full_path <- file.path(de_dir, f)
  resample_index <- gsub("pvals_", "", gsub(".csv", "", f, fixed = T), fixed = T)
  resample_index <- as.integer(resample_index)
  d <- read_csv(full_path)
  d <- mutate(d, resample_index = resample_index)
  return(d)
})

de_df <- bind_rows(dfs)
rm(dfs)

robust_df <- de_df %>%
  group_by(gene) %>%
  summarise(prop_sig = mean(q_val < alpha))

robust_genes <- robust_df %>%
  filter(prop_sig > (1 - alpha)) %>%
  extract2("gene")

load("data/sce_trapnell.Rdata")

sce <- sce[robust_genes, ]

save(sce, file = "data/resamples/sce_trapnell_robust.Rdata")

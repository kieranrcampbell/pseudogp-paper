library(scater)
library(readr)
library(dplyr)
library(matrixStats)
library(reshape2)
library(caret)

set.seed(123L)

theme_set(theme_bw())

alpha <- 0.05

resample_dir <- "data/resamples/diffexpr"

resample_files <- dir(resample_dir)

resample_qvals <- sapply(resample_files, function(f) {
  csv <- read_csv(file.path(resample_dir, f))
  csv$q_val
})

gene_names <- read_csv(file.path(resample_dir, resample_files[1]))$gene

gene_df <- data_frame(gene = gene_names,
                      p_sig = rowMeans(resample_qvals < alpha)) %>%
  mutate(is_sig = p_sig > (1 - alpha))

to_sample <- createDataPartition(gene_df$is_sig, p = 0.1)$Resample1

genes_to_use <- gene_names[to_sample]

load("data/resamples/sce_trapnell_resamples.Rdata")

sce <- sce[genes_to_use, ]

save(sce, file = "data/resamples/sce_trapnell_gplvm.Rdata")

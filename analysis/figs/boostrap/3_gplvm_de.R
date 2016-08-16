library(scater)
library(readr)
library(dplyr)
library(embeddr)

args <- commandArgs(trailingOnly = TRUE)
which_bootstrap <- as.integer(args[1])

output_file <- paste0("data/bootstrap/gplvm_de/de_", which_bootstrap, ".csv")

gplvm_pseudotimes <- read_csv("data/bootstrap/gplvm_pseudotimes.csv")

pst <- gplvm_pseudotimes[[which_bootstrap]]

load("data/sce_trapnell.Rdata")
sce@featureControlInfo <- AnnotatedDataFrame()

n_cells_exprs <- rowSums(exprs(sce) > sce@lowerDetectionLimit)
genes_to_use <- n_cells_exprs > (0.1 * ncol(sce)) # select genes expressed in at least 10% of cells
sce <- sce[genes_to_use,]

sce$pseudotime <- pst

de_test <- pseudotime_test(sce, n_cores = 1)

write_csv(de_test, output_file)


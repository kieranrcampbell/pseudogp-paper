library(scater)
library(readr)
library(dplyr)
library(embeddr)

args <- commandArgs(trailingOnly = TRUE)
which_bootstrap <- as.integer(args[1])

output_file <- paste0("data/bootstrap/bootstrapped_de/de_", which_bootstrap, ".csv")

which_cells <- read_csv("data/bootstrap/which_cells.csv")
bootstrapped_pseudotimes <- read_csv("data/bootstrap/bootstrapped_pseudotimes.csv")

cells <- which_cells[[which_bootstrap]]
pst <- bootstrapped_pseudotimes[[which_bootstrap]]

load("data/sce_trapnell.Rdata")
sce@featureControlInfo <- AnnotatedDataFrame()

n_cells_exprs <- rowSums(exprs(sce) > sce@lowerDetectionLimit)
genes_to_use <- n_cells_exprs > (0.1 * ncol(sce)) # select genes expressed in at least 10% of cells
sce <- sce[genes_to_use,]

new_exprs <- exprs(sce)[, cells]
colnames(new_exprs) <- paste0("Cell", 1:ncol(new_exprs))
new_pdata <- pData(sce)[cells,]
rownames(new_pdata) <- colnames(new_exprs)
new_sceset <- newSCESet(exprsData = new_exprs, phenoData = AnnotatedDataFrame(new_pdata))

new_sceset$pseudotime <- pst

de_test <- pseudotime_test(new_sceset, n_cores = 1)

write_csv(de_test, output_file)
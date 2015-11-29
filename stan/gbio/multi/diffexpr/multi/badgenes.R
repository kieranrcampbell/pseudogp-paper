## choose bad genes

library(readr)
library(dplyr)

x <- read_csv("multi/dfs.csv")

base_dir <- "~/mount/"
source(paste0(base_dir, "GP/pseudogp2/stan/diffexpr/monocle/prep_data.R"))
sce <- load_data(base_dir)
rownames(x) <- featureNames(sce)
x$gene <- fData(sce)$gene_short_name

ggplot(x) + geom_point(aes(x = psig, y = qval))

bad_genes <- filter(x, psig < 0.1, qval < 0.1)
bad_genes_sort <- arrange(bad_genes, qval)


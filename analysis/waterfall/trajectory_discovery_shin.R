#' # Trajectory discovery in the Shin et al. (2015) dataset
#' ### Kieran Campbell <kieran.campbell@sjc.ox.ac.uk>
#' 
#' To turn this into markdown, run `knitr::spin("trajectory_discovery_shin.R")` from
#' within R.

#' ### Loading the data
#+ setup
library(readxl)
library(scater)
library(matrixStats)
library(rhdf5)
library(embeddr)

base_dir <- "~/mount"
data_dir <- file.path(base_dir, "datasets/waterfall")
output_sce <- file.path(base_dir, "pseudogp-paper/data/sce_waterfall.Rdata")
h5outfile <- file.path(base_dir, "pseudogp-paper/data/waterfall_embeddings.h5")


#' We need to use `read_excel` from the excellent `readxl` package to get the data:
#' 
#+ readxl
fname <- dir(data_dir)[grep("mmc7", dir(data_dir))]
excel_file <- file.path(data_dir, fname)
x <- read_excel(excel_file)
x <- as.data.frame(x)

#' This requires a bunch of munging to get the gene expression matrix separate from the
#' pseudotime vector (as assigned by waterfall)
cellnames <- as.character(x[1,-1])
genenames <- x[-c(1,2, nrow(x)),1]
w_pseudotimes <- as.numeric(x[2,-1])
expr <- as.matrix(x[-c(1,2, nrow(x)),-1])
expr <- t(apply(expr, 1, as.numeric))

colnames(expr) <- cellnames
rownames(expr) <- genenames

#' Now turn into a `scater` object:
#+ to-scater
pd <- new("AnnotatedDataFrame", data.frame(wpst = w_pseudotimes))
rownames(pd) <- cellnames

sce <- newSCESet(expr = expr, phenoData = pd, logExprsOffset = 1)

ercc <- grep("ERCC", featureNames(sce), ignore.case = TRUE)
sce <- calculateQCMetrics(sce)


#' A lot of the calculations involve using top 195 most variable genes, which
#' we can easily calculate:
#+ most-var-genes
rv <- matrixStats::rowVars(exprs(sce))
ntop <- 195
feature_set <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]

#' ### Representations
#' 
#' To get a PCA representation, we use the top 195 genes but *don't* scale the features. This makes
#' it look most like 
#+ pca-plt, fig.width=6, fig.height=5
sce <- plotPCA(sce, colour_by = "wpst", ntop = ntop, return_SCESet = TRUE, scale_features = FALSE)
Xpca <- redDim(sce)

#' We can also get a tSNE representation. To make it look most 'trajectory-like', we use
#' a perplexity of 3:
#+ tsne-rep, fig.width=6, fig.height=5
set.seed(12)
sce <- plotTSNE(sce, colour_by = "wpst", ntop = ntop, return_SCESet = TRUE, perplexity = 3)
Xtsne <- redDim(sce)

#' Finally to get the laplacian eigenmaps embedding, we use `embeddr` with a euclidean
#' metric and 30 nearest neighbours, using only the top 195 most variable genes
#+ embeddr-plt, fig.width=6, fig.height=5
sce <- embeddr(sce, nn = 30, metric = "euclidean", genes_for_embedding = feature_set)
sce <- fit_pseudotime(sce)
plotReducedDim(sce, colour_by = "wpst")
Xle <- redDim(sce)
t_gt <- pseudotime(sce)

#' We can now save the `SCESet` and the three embeddings, as well as the principal curve 'ground truth'

#+ save_sceset
save(sce, file = output_sce)

h5createFile(h5outfile)
h5write(Xpca, h5outfile, "Xpca")
h5write(Xtsne, h5outfile, "Xtsne")
h5write(Xle, h5outfile, "Xle")
h5write(t_gt, h5outfile, "t_gt")
h5write(w_pseudotimes, h5outfile, "wpst")




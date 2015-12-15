library(readxl)
library(scater)
library(matrixStats)
library(rhdf5)
library(embeddr)

base_dir <- "~/mount"
data_dir <- file.path(base_dir, "datasets/waterfall")

fname <- dir(data_dir)[grep("mmc7", dir(data_dir))]
excel_file <- file.path(data_dir, fname)
x <- read_excel(excel_file)
x <- as.data.frame(x)
cellnames <- as.character(x[1,-1])
genenames <- x[-c(1,2, nrow(x)),1]
w_pseudotimes <- as.numeric(x[2,-1])
expr <- as.matrix(x[-c(1,2, nrow(x)),-1])
expr <- t(apply(expr, 1, as.numeric))

colnames(expr) <- cellnames
rownames(expr) <- genenames

pd <- new("AnnotatedDataFrame", data.frame(wpst = w_pseudotimes))
rownames(pd) <- cellnames

sce <- newSCESet(expr = expr, phenoData = pd, logExprsOffset = 1)

ercc <- grep("ERCC", featureNames(sce), ignore.case = TRUE)
sce <- calculateQCMetrics(sce)

rv <- matrixStats::rowVars(exprs(sce))
ntop <- 195
feature_set <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]

set.seed(12)
sce <- plotPCA(sce, colour_by = "wpst", ntop = ntop, return_SCESet = TRUE, scale_features = FALSE)
Xpca <- redDim(sce)
sce <- plotTSNE(sce, colour_by = "wpst", ntop = ntop, return_SCESet = TRUE, perplexity = 3)
Xtsne <- redDim(sce)

sce <- embeddr(sce, nn = 30, metric = "euclidean", genes_for_embedding = feature_set)
sce <- fit_pseudotime(sce)
Xle <- redDim(sce)
t_gt <- pseudotime(sce)

# save SCESet
save(sce, file = file.path(base_dir, "GP/pseudogp2/stan/waterfall/waterfall.Rdata"))

h5outfile <- file.path(base_dir, "GP/pseudogp2/data/waterfall_embeddings.h5")

h5createFile(h5outfile)
h5write(Xpca, h5outfile, "Xpca")
h5write(Xtsne, h5outfile, "Xtsne")
h5write(Xle, h5outfile, "Xle")
h5write(t_gt, h5outfile, "t_gt")





library(readxl)
library(scater)
library(embeddr)

setwd("/net/isi-scratch/kieran/datasets/nose")
fname <- dir()[grep("mmc7", dir())]

x <- read_excel(fname)
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

save(sce, file = "ear.Rdata")
load("ear.Rdata")

ercc <- grep("ERCC", featureNames(sce), ignore.case = TRUE)
sce <- calculateQCMetrics(sce)

sce <- plotPCA(sce, colour_by = "wpst", ntop = 195, return_SCESet = TRUE, scale_features = FALSE)



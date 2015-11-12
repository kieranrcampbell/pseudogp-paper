library(matrixStats)
library(scater)
library(embeddr)
library(readr)
library(fastICA)
library(rhdf5)


setwd("/net/isi-scratch/kieran/datasets/ear/")

fname <- dir()[grep("TPM_Matrix", dir())]

x <- as.data.frame(read_delim(fname, "\t", quote = ''))

xmarker <- x[1:2,]
xgene <- x[3:nrow(x),]

clean_rownames <- function(s) gsub('\"', '', s)
rownames(xgene) <- sapply(xgene[,1], clean_rownames) 
xgene <- xgene[,-1]

## first transformation: anything < 1TPM -> 0, anything > 1TPM -> log2
xgene[xgene < 1] <- 0
xgene[xgene > 0] <- log2(xgene[xgene > 0])

ercc <- grep("ERCC", rownames(xgene), ignore.case = TRUE)

is_ute <- grep("Ute_P1", colnames(xgene))

xute <- xgene[,is_ute]
is_gfp <- grep("GFP", colnames(xute))
is_tom <- grep("tdTom", colnames(xute))
is_neg <- grep("Neg", colnames(xute))

marker <- rep(NA, ncol(xute))
marker[is_gfp] <- "GFP"
marker[is_tom] <- "Tom"
marker[is_neg] <- "Neg"

pd <- data.frame(marker = marker)
rownames(pd) <- colnames(xute)
pd <- new('AnnotatedDataFrame', pd)

sce <- newSCESet(exprsData = xute, logged = TRUE,
                 phenoData = pd, is_exprsData = xute > 0)

n_cells_exprs <- rowSums(is_exprs(sce))
sce <- sce[n_cells_exprs > 2, ]
ccv <- rowVars(exprs(sce)) / rowMeans(exprs(sce))^2 
sce <- sce[ccv >= 0.5, ]
save(sce, file = "sce.Rdata")

load("sce.Rdata")
ercc <- grep("ERCC", featureNames(sce), ignore.case = TRUE)

sce <- calculateQCMetrics(sce, feature_controls = ercc)

plotPhenoData(sce)

plotPCA(sce, colour_by = "marker", ntop = 195, ncomponents = 2)
plotPCA(sce, colour_by = "marker", ntop = 1000, ncomponents = 3)


plotPCA(sce, colour_by = "depth", ntop = 1000, ncomponents = 2)
plotPCA(sce, colour_by = "pct_dropout", ntop = 1000, ncomponents = 3)
plotPCA(sce, colour_by = "coverage")


sc <- embeddr(sce, nn = 20, metric = "cosine")
sc <- fit_pseudotime(sc)
plot_embedding(sc, color_by = "marker")
sc <- cluster_embedding(sc)
plot_embedding(sc, color_by = "marker")

pub_genes <- c("$Actb^", "Pou4f3", "Atoh1", "$Espn^", "Dkk3",
               "Myo7a", "Hes6", "Xirp2", "Hes1", "Gfi1", "Jag2", "Fscn2")
pub_inds <- grep(paste(pub_genes, collapse="|"), featureNames(sce), ignore.case = TRUE)

plot_in_pseudotime(sc[pub_inds,], color_by = "marker")
plot_pseudotime_model(sc[pub_inds,], color_by = "marker", line_color = "black")

xmarker <- xmarker[,-1]
ute_marker <- xmarker[,is_ute]
ute_marker <- data.frame(t(ute_marker))
names(ute_marker) <- c("GFP", "Tom")
ute_marker$pseudotime <- sc$pseudotime
ute_marker$marker <- marker
ggplot(ute_marker) + geom_point(aes(x = GFP, y = Tom, color = marker), size = 3)

um <- reshape2::melt(ute_marker, id.vars = "pseudotime", variable.name = "marker")
ggplot(um) + geom_point(aes(x = pseudotime, y = value)) + facet_wrap(~marker, scales = "free_y")

ica <- fastICA(exprs(sce), n.comp = 2)
redDim(sce) <- t(ica$A)

plotReducedDim(sce, colour_by = "marker")


# “cluster” cells as in paper ---------------------------------------------

sce <- plotPCA(sce, ntop = 195, return_SCESet = TRUE)

set.seed(1234)
km <- kmeans(redDim(sce), 6)
pData(sce)$cluster <- km$cluster
plotReducedDim(sce, colour_by = "cluster", shape_by = "marker")

trajectory_cells <- sce$cluster %in% c(5, 3, 6)
sct <- sce[, trajectory_cells]

pub_genes <- c("$Actb^", "Pou4f3", "Atoh1", "$Espn^", "Dkk3",
               "Myo7a", "Hes6", "Xirp2", "Hes1", "Gfi1", "Jag2", "Fscn2")
pub_inds <- grep(paste(pub_genes, collapse="|"), featureNames(sct), ignore.case = TRUE)

sct <- embeddr(sct, nn = 30, metric = "cosine")
sct <- fit_pseudotime(sct)
sct <- reverse_pseudotime(sct)
plot_embedding(sct, color_by = "marker")
plot_in_pseudotime(sct[pub_inds,], color_by = "marker")

save(sct, file = "sce_pst.Rdata")

h5file <- "/net/isi-scratch/kieran/GP/pseudogp2/data/ear_embeddings.h5"
h5createFile(h5file)
h5write(redDim(sct), file = h5file, "X")
h5write(pseudotime(sct), file = h5file, "t_gt")


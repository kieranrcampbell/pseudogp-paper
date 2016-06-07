
#' # Trajectory discovery in the Burns et al. (2015) dataset
#' ### Kieran Campbell <kieran.campbell@sjc.ox.ac.uk>

#' 
#' ## Introduction
#' The Burns et al. dataset contains single-cell RNA-seq data for 301 cells from the
#' utricilar and cochlear sensory epithelia of newborn mice. This document goes through
#' reading in the data and attempting to find a trajectory, with some discussion of why
#' clustering may be more appropriate.
#' 
#' To turn this into markdown, run `knitr::spin("burns.R")` from
#' within R.
#' 
#' ## Obtaining the data
#+ setup
library(matrixStats)
library(scater)
library(embeddr)
library(readr)
library(rhdf5)
library(reshape2)
library(cowplot)
library(rstan)
library(MCMCglmm)
library(coda)
library(ggplot2)
library(moments)
library(pseudogp)

pdf("figs/diagnostic/burns.pdf", width = 9)


#' For the purposes of this we've downloaded the TPM matrix, which can then
#' be easily read in using `read_delim` from the `readr` package:
#' 
#+ read-data, message = FALSE

www <- "http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE71982&format=file&file=GSE71982%5FRSEM%5FTPM%5FMatrix%2Etxt%2Egz"
filename <- "data/burns_raw.txt.gz"

if(!file.exists(filename)) {
  download.file(www, filename)
}

output_hdf5 <- "data/burns_embeddings.h5"
output_sce <- "data/sce_burns.Rdata"
pst_output_hdf5 <- "data/burns_pseudotime_traces.h5"

pdf("figs/diagnostic/burns.pdf")
x <- as.data.frame(read_delim(filename, "\t", quote = ''))

#' The first two rows of `x` (as far as we can tell) contain fluorescent intensity information about
#' EGFP and tdTom so we can assign those to a separate data frame and all the transcription data
#' to a separate one, and tidy up the rownames of the genes

#+ sep-matrices
clean_rownames <- function(s) gsub('\"', '', s)
rownames(x) <- sapply(x[,1], clean_rownames) 
x <- x[,-1]

xmarker <- x[1:2,]
xgene <- x[3:nrow(x),]

#' The methods part of the section perform rather bizarre gene normalisation: 
#' > For all analyses, nTPMs < 1 were set to zero, and TPM > 1 were transformed to log2(nTPM) space.
#' 
#' The authors define nTPM as TPM that's been scaled by its geometric mean. The DESeq(2) papers explain
#' quite clearly you shouldn't do this as TPM are already normalised. For the sake of brevity we'll
#' apply the log normalisation and leave it there:

#+ log-normalise
xgene[xgene < 1] <- 0
xgene[xgene > 0] <- log2(xgene[xgene > 0])

#' The paper mainly looks at P1 utricular cells, which we'll assume are denoted by "Ute_P1" 
#' in the cell names, so we'll use those for now
#+ ute-only
is_ute <- grep("Ute_P1", colnames(xgene))
xute <- xgene[,is_ute]

#' This gives us 160 cells rather than the 158 claimed in the paper, but since none 
#' have a differently named format from the rest it's impossible to tell which the two
#' outliers are, so we'll keep all 160 for now.

#' Let's have a look at what these markers look like
#+ markers, fig.width = 7, fig.height = 5
is_gfp <- grep("GFP", colnames(xute))
is_tom <- grep("tdTom", colnames(xute))
is_neg <- grep("Neg", colnames(xute))

marker <- rep(NA, ncol(xute))
marker[is_gfp] <- "GFP"
marker[is_tom] <- "Tom"
marker[is_neg] <- "Neg"

xute_marker <- xmarker[,is_ute]
df_marker <- data.frame(t(xute_marker), designated_marker = marker)
dfm <- melt(df_marker, id.vars = "designated_marker", variable.name = "measured_marker")
ggplot(dfm) + geom_boxplot(aes(color = measured_marker, y = value, x = designated_marker)) + 
  scale_y_log10()

#' So this roughly checks out, though there's clearly a lot of coexpression of EGFP and tdTom
#' in the cells designated as tdTom. Time to turn it into a `Scater` object:

#+ to-scater
pd <- data.frame(marker = marker)
rownames(pd) <- colnames(xute)
pd <- new('AnnotatedDataFrame', pd)

sce <- newSCESet(exprsData = xute, logged = TRUE,
                 phenoData = pd, is_exprsData = xute > 0)
ercc <- grep("ERCC", featureNames(sce), ignore.case = TRUE)
sce <- calculateQCMetrics(sce, feature_controls = ercc)

#' ## Identifying trajectories

#' We can make a basic PCA plot using the 'top' 195 genes (in scater this is defined as most variable).
#' We set `r scale_features = FALSE` as this is clearly what's done in the paper.
#+ first-pca, fig.width=6, fig.height = 4
sce <- plotPCA(sce, colour_by = "marker", ntop = 195, 
               ncomponents = 2, return_SCESet = TRUE, scale_features = FALSE)
redDim(sce) <- redDim(sce)[,1:2]

#' Now this looks like figure 4a flipped across the x-axis. Does this make sense with
#' respect to fluorescent markers? The paper states
#' > all HCs in both the utricle and cochlea express tdTomato driven by Gfi1Cre, whereas 1% of SCs are
#' > tdTomato+ve. In addition, most SCs in the utricle and cochlea express high levels of green fluorescent protein (GFP), 
#' > driven by the Lfng promoter. Many utricular HCs and some cochlear HCs also express GFP, but generally at lower levels. 
#' > Striolar SCs and transitional epithelial cells (TECs) located at the border between sensory and non-sensory regions in the utricle, 
#' > and inner pillar cells and non-sensory cells (NSCs) in the cochlea express low or undetectable levels of both fluorescent proteins.
#' 
#' This would imply the blue GFP+ve cluster are SCs, the green tdTomato cells are HCs and the orange neg cells
#' are TECs, consistent with the notion that this is fig 4a mirrored. To remove any doubt, let's plot
#' the three marker genes in figure 4c:
#+ plot-markers, fig.width=10, fig.height=3
Lfng_plt <- plotReducedDim(sce, colour_by = "marker", size_by = "Lfng")
Cdh4_plt <- plotReducedDim(sce, colour_by = "marker", size_by = "Cdh4")
Sox2_plt <- plotReducedDim(sce, colour_by = "marker", size_by = "Sox2")
plot_grid(Lfng_plt, Cdh4_plt, Sox2_plt, nrow = 1)

#' If you stare at these marker gene plots and refer to figure 4c it becomes evident that our PCA
#' is indeed that of figure 4a mirrored. We can therefore use k-means to cluster the cells into
#' 6 types and retain only those that correspond to the green -> yellow -> red SC to HC transition
#' in the original paper
#+ find-trajectory, fig.width=6, fig.height=4
set.seed(1234)
km <- kmeans(redDim(sce), 6)
pData(sce)$cluster <- km$cluster
plotReducedDim(sce, colour_by = "cluster", shape_by = "marker")

trajectory_cells <- sce$cluster %in% c(3,5,6)
sct <- sce[, trajectory_cells]


#' Let's have a go at finding a trajectory with `embeddr`, using the top 195 most
#' variable genes with the standard correlation metric and 20 nearest neighbours:
#+ quick-traj, fig.width=5, fig.height=5
rv <- matrixStats::rowVars(exprs(sct))
ntop <- 195
feature_set <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
plot_embedding(embeddr(sct, nn = 25, genes_for_embedding = feature_set))

#' There seems to be quite distinct clustering between the two datasets (laplacian eigenmap
#' embeddings that look like this imply the graph almost has two connected components and we
#' might as well use spectral clustering). What about using the entire gene set?
#+ embeddr-entire, fig.width=5, fig.height=5
plot_embedding(embeddr(sct, nn = 25))

#' There's now a more clearly defined trajectory. We can try and fit a principal curve to it:
#+ fit-curve, fig.width = 5, fig.height = 5
sct <- embeddr(sct, nn = 25)
sct <- sct[,redDim(sct)[,2] < 0.24] # remove outliers
sct <- fit_pseudotime(sct)
Xle <- redDim(sct)[,1:2] # Laplacian eigenmaps embedding we'll use
plot_embedding(sct)


#' Let's do a final sanity check and make sure the direction of the genes is the same as that
#' in figure 5b (or at least reversed, remember pseudotime is pseudo!).

#+ final-sanity, fig.width = 9, fig.height=7
pub_genes <- c("$Actb^", "Pou4f3", "Atoh1", "$Espn^", "Dkk3",
               "Myo7a", "Hes6", "Xirp2", "Hes1", "Gfi1", "Jag2", "Fscn2")
pub_inds <- grep(paste(pub_genes, collapse="|"), featureNames(sct), ignore.case = TRUE)
plot_in_pseudotime(reverse_pseudotime(sct[pub_inds,]), color_by = "marker")

#' Looking good! But is it actually a trajectory, or is it just an artefact of us forcing
#' a k-nearest-neighbour graph onto it?
#' 
#' Let's plot the PCA representation again:
#+ pca-2, fig.width=6, fig.height=6
sct <- plotPCA(sct, colour_by = "pseudotime", scale_features = FALSE, return_SCESet = TRUE)
Xpca <- redDim(sct)[, 1:2]

#' Again, there's no *real* trajectory structure here. What about t-SNE? Let's be fair and plot over 
#' a range of perplexities in case a certain one gives us something more amenable to 
#' trajectory discovery:
#+ tsnes, fig.width=9, fig.height=10
set.seed(123) # this is tsne so we need a seed
perplexities <- c(1,2,5,10,20)
tsne_plts <- lapply(perplexities, function(p) plotTSNE(sct, colour_by = "marker", scale_features = FALSE, perplexity = p))
plot_grid(plotlist = tsne_plts, nrow = 3, labels = as.character(perplexities))

#' At very low perplexity we might be able to fit a trajectory, but really these cells cluster rather 
#' beautifully on fluorescent marker expression.
#' 
#' Now we want to save the embeddings for use in the rest of the analysis. We'll choose
#' a perplexity of 2 for the tSNE plotting as this appears to give the most 
#' "trajectory-like" structure.
#+ get-tsne-x
set.seed(123)
sct <- plotTSNE(sct, colour_by = "pseudotime", shape_by = "marker", 
                scale_features = FALSE, return_SCESet = TRUE, perplexity = 2)
Xtsne <- redDim(sct)

#' Now we can save everything to HDF5. We'll call the pseudotime fit
#' from `embeddr` as `t_gt`:
#+ save-all
t_gt <- pseudotime(sct)
h5createFile(output_hdf5)
h5write(Xle, output_hdf5, "Xle")
h5write(Xpca, output_hdf5,"Xpca")
h5write(Xtsne, output_hdf5, "Xtsne")
h5write(t_gt, output_hdf5, "t_gt")

#' and save the `SCESet` object too:
#+ save-sce
sce <- sct
save(sce, file = output_sce)

#+ Pseudotime fitting

X <- Xle
#X <- apply(X, 2, function(x) (x - mean(x)) / sd(x))

#' and quickly plot to make sure it looks right
#+ quick-plot, fig.width=5, fig.height = 5
ggplot(data.frame(X, t_gt)) + 
  geom_point(aes(x = component_1, y = component_2, color = t_gt)) + theme_bw()

#' ### Fit the pseudotime
#+ pseudo-fit
fit <- fitPseudotime(X, initialise_from = "pca", 
                     smoothing_alpha = 14, smoothing_beta = 1, seed = 123,
                     iter = 12000, thin = 6)

#' Diagnostics
#+ plot-diagnostic
plotDiagnostic(fit)


#' Plot posterior mean curves
#+ posteriorcplt, fig.width=6, fig.height=5
posteriorCurvePlot(X, fit, posterior_mean = TRUE)

#' Posterior boxplot of pseudotime distribution
#+ posteriorbplt, fig.width=7, fig.height=5
posteriorBoxplot(fit)

#' Because we're dealing with stan objects we can very easily extract the posterior samples:
#+ extract samples
pst <- extract(fit, "t")
tmcmc <- mcmc(pst$t)

smcmc <- mcmc(extract(fit, "sigma")$sigma[,1,]) # need to slice middle index to get single rep.
lmcmc <- mcmc(extract(fit, "lambda")$lambda[,1,])

#' ...and save them to HDF5
#+ save-hdf5
if(!file.exists(pst_output_hdf5)) h5createFile(pst_output_hdf5)
h5write(X, pst_output_hdf5, "X")
h5write(pst$t, pst_output_hdf5, "pst")
h5write(as.matrix(smcmc), pst_output_hdf5, "sigma")
h5write(as.matrix(lmcmc), pst_output_hdf5, "lambda")

dev.off()


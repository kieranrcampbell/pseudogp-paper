#' ## Trajectory discovery for Trapnell et al (2014)
#' ### Kieran Campbell <kieran.campbell@sjc.ox.ac.uk>
#' 
#' To turn this into markdown, run `knitr::spin("trajectory_discovery_trapnell.R")` from
#' within R.

#' Here we show how `embeddr` (= spectral embedding + principal curves) can be used for pseudotemporal ordering of single-cell gene expression data using the [monocle](http://cole-trapnell-lab.github.io/monocle-release/) dataset. This uses the `HSMMSingleCell` dataset that is bundled with monocle.

#+ load-all, message=FALSE, warning=FALSE
library(monocle) ## for monocle data
library(devtools) ## for package development
library(ggplot2)
library(rhdf5)
library(scater) 
library(embeddr)
library(readr)
library(plyr)

base_dir <- "~/mount"

output_hdf5 <- file.path(base_dir, "pseudogp-paper/data/trapnell_embeddings.h5")
output_sce <- file.path(base_dir, "pseudogp-paper/data/sce_monocle.Rdata")

#' First we create the `SCESet` using the data from the `HSMMSingleCell` package.
#' This is a bit fiddly since HSMMSingleCell changed format recently
#+ create-sce, message=FALSE, warning=FALSE
sce <- NULL
hsmm_data_available <- data(package='HSMMSingleCell')$results[,3]
if("HSMM" %in% hsmm_data_available) {
  data(HSMM)
  sce <- fromCellDataSet(HSMM, exprs_values = 'fpkm')
} else if("HSMM_expr_matrix" %in% hsmm_data_available) {
  data(HSMM_expr_matrix)
  data(HSMM_gene_annotation)
  data(HSMM_sample_sheet)

  pd <- new('AnnotatedDataFrame', data = HSMM_sample_sheet)
  fd <- new('AnnotatedDataFrame', data = HSMM_gene_annotation)
  sce <- newSCESet(fpkmData = HSMM_expr_matrix, phenoData = pd, featureData = fd)
} else {
  stop('No recognised data types in HSMMSingleCell')
}

## add cell_id to HSMM to play nicely with dplyr
sce$cell_id <- rownames(pData(sce))

#' This dataset seems to work best in log10 so we'll rescale to use that:
#+ rescale
exprs(sce) <- log10(fpkm(sce) + 1)

#' ### Selecting genes for the embedding
#' We'll fit a CV2 - mean curve to pick out highly variable genes for the embedding

#+ cv2, fig.width=6, fig.height=5
x <- t(exprs(sce)) # t(log10(exprs(sce) + 1))
x_mean <- colMeans(x)
x_var <- apply(x, 2, var)
genes_for_fit <- x_mean > 0.3
CV2 <- x_var[genes_for_fit] / (x_mean[genes_for_fit])^2
df_fit <- data.frame(m = x_mean[genes_for_fit], CV2 = CV2)
fit_loglin <- nls(CV2 ~ a * 10^(-k * m), data = df_fit, start=c(a=5, k=1)) 
ak <- coef(fit_loglin)
f <- function(x) ak[1] * 10^(-ak[2] * x) 
genes_for_embedding <- (CV2 > 4 * predict(fit_loglin))
df_fit$for_embedding <- as.factor(genes_for_embedding)
df_fit$gene_names <- names(genes_for_embedding)

ggplot(df_fit, aes(x=m, y=CV2, color = for_embedding)) + geom_point() +
  theme_bw() + xlab('Mean') + ylab('CV2') + 
  stat_function(fun=f, color='black')

#' Next we take the log10 of the dataset (using a pseudocount of 1) and fit the embedding using the `embeddr` function using the default settings:

#+ embedding_highmag, fig.width=7.5, fig.height=4.5
set.seed(123)
#sce <- fromCellDataSet(HSMM, use_exprs_as = "fpkm")

gene_indices <- match(names(which(genes_for_embedding)), featureNames(sce))
sce <- embeddr(sce, genes_for_embedding = gene_indices)

pData(sce)$long_state <- plyr::mapvalues(pData(sce)$State, from=1:3,
                                            to=c('Proliferating cell',
                                                 'Differentiating myoblast',
                                                 'Interstitial mesenchymal cell'))

plot_embedding(sce, color_by = 'long_state')



#' We can also cluster the embedding using kmeans to get rid of contaminating cells
#+ clust_emb, fig.width=7.5, fig.height=4.5
sce <- cluster_embedding(sce, k = 3)

sce_tmp <- sce
phenoData(sce_tmp)$cluster <- plyr::mapvalues(pData(sce_tmp)$cluster, from=c(3, 1, 2),
                                            to=c(1,2,3))
phenoData(sce_tmp)$cluster <- plyr::mapvalues(pData(sce_tmp)$cluster, from=1:3,
                                            to=c('Interstitial mesenchymal cell',
                                                 'Proliferating cell',
                                                 'Differentiating myoblast'))

plot_embedding(sce_tmp)


#' ### Pseudotime fitting
#' In  Trapnell et al. paper they show that groups 1 & 3 correspond 
#' to differentiating cells while group 2 is contamination. 
#' We can separate off groups 1 & 3, fit the pseudotime trajectories and plot:

#+ fit-pseudotime
sce_23 <- sce[, pData(sce)$cluster %in% c(1,2)]
## remove a couple of outlier cells
sce_23 <- sce_23[,redDim(sce_23)[,1] > 0]
sce_23 <- fit_pseudotime(sce_23)
plot_embedding(sce_23)
Xle <- redDim(sce_23)

#' Let's also look at PCA and tSNE
#+ pca, fig.width=6, fig.height=5
sce_23 <- plotPCA(sce_23, colour_by = "pseudotime", return_SCESet = TRUE)
Xpca <- redDim(sce_23)

#' And tSNE:
#+ tsne, fig.width=6, fig.height=5
set.seed(1234)
sce_23 <- plotTSNE(sce_23, colour_by = "pseudotime", perplexity = 3, return_SCESet = TRUE)
Xtsne <- redDim(sce_23)

#' Then save to HDF5 in the usual way:
#' 
#+ to-hdf5
h5createFile(output_hdf5)
h5write(Xle, output_hdf5, "Xle")
h5write(pseudotime(sce_23), output_hdf5, "t_gt")
h5write(Xpca, output_hdf5, "Xpca")
h5write(Xtsne, output_hdf5, "Xtsne")

save(sce_23, file = output_sce)
  
  

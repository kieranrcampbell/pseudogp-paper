# Trajectory discovery in the Burns et al. (2015) dataset
### Kieran Campbell <kieran.campbell@sjc.ox.ac.uk>

## Introduction
The Burns et al. dataset contains single-cell RNA-seq data for 301 cells from the
utricilar and cochlear sensory epithelia of newborn mice. This document goes through
reading in the data and attempting to find a trajectory, with some discussion of why
clustering may be more appropriate.

To turn this into markdown, run `knitr::spin("trajectory_discovery_burns.R")` from
within R.

## Obtaining the data


```r
library(matrixStats)
library(scater)
library(embeddr)
library(readr)
library(rhdf5)
library(reshape2)
library(cowplot)
```

For the purposes of this we've downloaded the TPM matrix, which can then
be easily read in using `read_delim` from the `readr` package:



```r
base_dir <- "~/mount"
data_dir <- file.path(base_dir, "datasets/ear/")
output_hdf5 <- file.path(base_dir, "pseudogp-paper/data/ear_embeddings.h5")
output_sce <- file.path(base_dir, "pseudogp-paper/data/sce_ear.Rdata")

fname <- dir(data_dir)[grep("TPM_Matrix", dir(data_dir))]
data_file <- file.path(data_dir, fname)

x <- as.data.frame(read_delim(data_file, "\t", quote = ''))
```

The first two rows of `x` (as far as we can tell) contain fluorescent intensity information about
EGFP and tdTom so we can assign those to a separate data frame and all the transcription data
to a separate one, and tidy up the rownames of the genes


```r
clean_rownames <- function(s) gsub('\"', '', s)
rownames(x) <- sapply(x[,1], clean_rownames) 
x <- x[,-1]

xmarker <- x[1:2,]
xgene <- x[3:nrow(x),]
```

The methods part of the section perform rather bizarre gene normalisation: 
> For all analyses, nTPMs < 1 were set to zero, and TPM > 1 were transformed to log2(nTPM) space.

The authors define nTPM as TPM that's been scaled by its geometric mean. The DESeq(2) papers explain
quite clearly you shouldn't do this as TPM are already normalised. For the sake of brevity we'll
apply the log normalisation and leave it there:


```r
xgene[xgene < 1] <- 0
xgene[xgene > 0] <- log2(xgene[xgene > 0])
```

The paper mainly looks at P1 utricular cells, which we'll assume are denoted by "Ute_P1" 
in the cell names, so we'll use those for now


```r
is_ute <- grep("Ute_P1", colnames(xgene))
xute <- xgene[,is_ute]
```

This gives us 160 cells rather than the 158 claimed in the paper, but since none 
have a differently named format from the rest it's impossible to tell which the two
outliers are, so we'll keep all 160 for now.
Let's have a look at what these markers look like


```r
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
```

```
## Warning: Removed 77 rows containing non-finite values (stat_boxplot).
```

![plot of chunk markers](figure/markers-1.png) 

So this roughly checks out, though there's clearly a lot of coexpression of EGFP and tdTom
in the cells designated as tdTom. Time to turn it into a `Scater` object:


```r
pd <- data.frame(marker = marker)
rownames(pd) <- colnames(xute)
pd <- new('AnnotatedDataFrame', pd)

sce <- newSCESet(exprsData = xute, logged = TRUE,
                 phenoData = pd, is_exprsData = xute > 0)
ercc <- grep("ERCC", featureNames(sce), ignore.case = TRUE)
sce <- calculateQCMetrics(sce, feature_controls = ercc)
```

## Identifying trajectories
We can make a basic PCA plot using the 'top' 195 genes (in scater this is defined as most variable).
We set  as this is clearly what's done in the paper.


```r
sce <- plotPCA(sce, colour_by = "marker", ntop = 195, 
               ncomponents = 2, return_SCESet = TRUE, scale_features = FALSE)
```

![plot of chunk first-pca](figure/first-pca-1.png) 

Now this looks like figure 4a flipped across the x-axis. Does this make sense with
respect to fluorescent markers? The paper states
> all HCs in both the utricle and cochlea express tdTomato driven by Gfi1Cre, whereas 1% of SCs are
> tdTomato+ve. In addition, most SCs in the utricle and cochlea express high levels of green fluorescent protein (GFP), 
> driven by the Lfng promoter. Many utricular HCs and some cochlear HCs also express GFP, but generally at lower levels. 
> Striolar SCs and transitional epithelial cells (TECs) located at the border between sensory and non-sensory regions in the utricle, 
> and inner pillar cells and non-sensory cells (NSCs) in the cochlea express low or undetectable levels of both fluorescent proteins.

This would imply the blue GFP+ve cluster are SCs, the green tdTomato cells are HCs and the orange neg cells
are TECs, consistent with the notion that this is fig 4a mirrored. To remove any doubt, let's plot
the three marker genes in figure 4c:


```r
Lfng_plt <- plotReducedDim(sce, colour_by = "marker", size_by = "Lfng")
Cdh4_plt <- plotReducedDim(sce, colour_by = "marker", size_by = "Cdh4")
Sox2_plt <- plotReducedDim(sce, colour_by = "marker", size_by = "Sox2")
plot_grid(Lfng_plt, Cdh4_plt, Sox2_plt, nrow = 1)
```

![plot of chunk plot-markers](figure/plot-markers-1.png) 

If you stare at these marker gene plots and refer to figure 4c it becomes evident that our PCA
is indeed that of figure 4a mirrored. We can therefore use k-means to cluster the cells into
6 types and retain only those that correspond to the green -> yellow -> red SC to HC transition
in the original paper


```r
set.seed(1234)
km <- kmeans(redDim(sce), 6)
pData(sce)$cluster <- km$cluster
plotReducedDim(sce, colour_by = "cluster", shape_by = "marker")
```

![plot of chunk find-trajectory](figure/find-trajectory-1.png) 

```r
trajectory_cells <- sce$cluster %in% c(3,5,6)
sct <- sce[, trajectory_cells]
```

Let's have a go at finding a trajectory with `embeddr`, using the top 195 most
variable genes with the standard correlation metric and 20 nearest neighbours:


```r
rv <- matrixStats::rowVars(exprs(sct))
ntop <- 195
feature_set <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
plot_embedding(embeddr(sct, nn = 25, genes_for_embedding = feature_set))
```

![plot of chunk quick-traj](figure/quick-traj-1.png) 

There seems to be quite distinct clustering between the two datasets (laplacian eigenmap
embeddings that look like this imply the graph almost has two connected components and we
might as well use spectral clustering). What about using the entire gene set?


```r
plot_embedding(embeddr(sct, nn = 25))
```

![plot of chunk embeddr-entire](figure/embeddr-entire-1.png) 

There's now a more clearly defined trajectory. We can try and fit a principal curve to it:


```r
sct <- embeddr(sct, nn = 25)
sct <- sct[,redDim(sct)[,2] < 0.3] # cheekily remove an outlier
sct <- fit_pseudotime(sct)
Xle <- redDim(sct) # Laplacian eigenmaps embedding we'll use
plot_embedding(sct)
```

![plot of chunk fit-curve](figure/fit-curve-1.png) 

Let's do a final sanity check and make sure the direction of the genes is the same as that
in figure 5b (or at least reversed, remember pseudotime is pseudo!).


```r
pub_genes <- c("$Actb^", "Pou4f3", "Atoh1", "$Espn^", "Dkk3",
               "Myo7a", "Hes6", "Xirp2", "Hes1", "Gfi1", "Jag2", "Fscn2")
pub_inds <- grep(paste(pub_genes, collapse="|"), featureNames(sct), ignore.case = TRUE)
plot_in_pseudotime(reverse_pseudotime(sct[pub_inds,]), color_by = "marker")
```

![plot of chunk final-sanity](figure/final-sanity-1.png) 

Looking good! But is it actually a trajectory, or is it just an artefact of us forcing
a k-nearest-neighbour graph onto it?

Let's plot the PCA representation again:


```r
sct <- plotPCA(sct, colour_by = "pseudotime", scale_features = FALSE, return_SCESet = TRUE)
```

![plot of chunk pca-2](figure/pca-2-1.png) 

```r
Xpca <- redDim(sct)
```

Again, there's no *real* trajectory structure here. What about t-SNE? Let's be fair and plot over 
a range of perplexities in case a certain one gives us something more amenable to 
trajectory discovery:


```r
set.seed(123) # this is tsne so we need a seed
perplexities <- c(1,2,5,10,20,30)
tsne_plts <- lapply(perplexities, function(p) plotTSNE(sct, colour_by = "marker", scale_features = FALSE, perplexity = p))
plot_grid(plotlist = tsne_plts, nrow = 3, labels = as.character(perplexities))
```

![plot of chunk tsnes](figure/tsnes-1.png) 

At very low perplexity we might be able to fit a trajectory, but really these cells cluster rather 
beautifully on fluorescent marker expression.

Now we want to save the embeddings for use in the rest of the analysis. We'll choose
a perplexity of 2 for the tSNE plotting as this appears to give the most 
"trajectory-like" structure.


```r
set.seed(123)
sct <- plotTSNE(sct, colour_by = "pseudotime", shape_by = "marker", 
                scale_features = FALSE, return_SCESet = TRUE, perplexity = 2)
```

![plot of chunk get-tsne-x](figure/get-tsne-x-1.png) 

```r
Xtsne <- redDim(sct)
```

Now we can save everything to HDF5. We'll call the pseudotime fit
from `embeddr` as `t_gt`:


```r
h5createFile(output_hdf5)
```

```
## [1] TRUE
```

```r
h5write(Xle, output_hdf5, "Xle")
h5write(Xpca, output_hdf5,"Xpca")
h5write(Xtsne, output_hdf5, "Xtsne")
h5write(pseudotime(sct), output_hdf5, "t_gt")
```

and save the `SCESet` object too:


```r
save(sct, file = output_sce)
```


library(monocle)
library(scater)
library(rhdf5)
library(readr)
library(coda)
library(MCMCglmm)
library(cowplot)
#library(pseudogp)
library(GGally)

# the pseudogp version does matrix standardizing
standardize <- function(t) (t - min(t))/(max(t) - min(t))


base_dir <- "~/mount/"
trace_file <- file.path(base_dir, "pseudogp-paper/data/monocle_stan_traces.h5")
sce_file <- file.path(base_dir, "pseudogp-paper/data/sce_monocle.Rdata")
embeddings_file <- file.path(base_dir, "pseudogp-paper/data/trapnell_embeddings.h5")
fig_path <- file.path(base_dir, "pseudogp-paper/analysis/figs/mstcomp/")

X <- h5read(embeddings_file, "Xle")

load(sce_file)
sce <- sce_23

n_cells_exprs <- rowSums(exprs(sce) > sce@lowerDetectionLimit)
genes_to_use <- n_cells_exprs > (0.1 * ncol(sce)) # select genes expressed in at least 10% of cells
sce <- sce[genes_to_use,]

# comparison of pseudotimes -----------------------------------------------

redDim(sce) <- X
#plotReducedDim(sce)
cds <- toCellDataSet(sce)
cds <- orderCells(cds)
cds_le_pseudotime <- standardize(cds$Pseudotime)

rv <- matrixStats::rowVars(exprs(sce))
ntop <- 500
feature_set <- 
  order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
cds <- setOrderingFilter(cds, featureNames(sce)[feature_set])
cds <- reduceDimension(cds, use_irlba = FALSE)
cds <- orderCells(cds)
#plot_spanning_tree(cds, color_by = "pseudotime")
monocle_pst <- cds$Pseudotime

tmcmc <- mcmc(h5read(trace_file, "pst"))
tmap <- posterior.mode(tmcmc)
intervals <- HPDinterval(tmcmc, prob = 0.95)

dft <- data.frame(monocle = monocle_pst, tmap = tmap, intervals, monocle_le = cds_le_pseudotime)

R2 <- signif(cor(dft$monocle, dft$tmap)^2, digits = 2)
R2text <- as.expression(substitute(italic(R)^2 == r, list(r = R2)))

cplt <- ggplot(dft) + geom_point(aes(x = monocle_pst, y = tmap)) + 
  geom_errorbar(aes(x = monocle_pst, ymin = lower, ymax = upper), alpha = 0.6, width = 1) +
  xlab("Monocle estimate") + ylab("GPLVM MAP estimate") +
  geom_text(x = 2, y = 0.9, label = as.character(R2text), parse = TRUE)

R2le <- signif(cor(dft$monocle_le, dft$tmap)^2, digits = 2)
R2textle <- as.expression(substitute(italic(R)^2 == r, list(r = R2le)))

leplt <- ggplot(dft) + geom_point(aes(x = monocle_le, y = tmap)) + 
  geom_errorbar(aes(x = monocle_le, ymin = lower, ymax = upper), alpha = 0.6, width = 0.05) +
  xlab("Monocle estimate") + ylab("GPLVM MAP estimate") +
  geom_text(x = 0.1, y = 0.9, label = as.character(R2textle), parse = TRUE)

corplts <- plot_grid(cplt, leplt, labels = c("ICA", "LE"), nrow = 1)

ggsave(file.path(fig_path, "S1_monocle_comparison.png"), corplts, width = 7, height = 3)


set.seed(123)
nsamples <- 30
ncells <- ncol(sce)
pseudosamples <- vector("list", ncells)
gt <- orderCells(cds)$Pseudotime
gt <- standardize(gt)

# did you ever see such awful R in your life
for(i in 1:nsamples) {
  chosen <- sample(ncells, round(0.8*ncells), replace=FALSE)
  cds_copy <- cds[,chosen]
  cds_copy@reducedDimS <- cds@reducedDimS[,chosen]
  cds_copy <- orderCells(cds_copy)
  pseudotimes <- standardize(cds_copy$Pseudotime)
  ## need to work out if we need to reverse the pseudotimes
  if(cor(gt[chosen], pseudotimes) < 0) pseudotimes <- 1 - pseudotimes
  
  for(j in 1:length(chosen)) {
    pseudosamples[[ chosen[j] ]] <- c(pseudosamples[[ chosen[j] ]], pseudotimes[j])
  }
}

pst_means <- sapply(pseudosamples, mean)
pst_sd <- sapply(pseudosamples, sd)
pst_95 <- 2*pst_sd


dplt <- ggplot(data.frame(pst_95)) + geom_density(aes(x = pst_95), fill = "black") +
  xlab(expression(2 * sigma)) + ylab("Density")

## boxplot
xp <- data.frame(t = do.call("c", pseudosamples))
nsamples <- sapply(pseudosamples, length)
cell <- lapply(1:length(nsamples), function(i) rep(i, nsamples[i]))
cell <- do.call("c", cell)
xp$cell <- as.factor(cell)

mbplt <- ggplot(xp) + geom_boxplot(aes(x = cell, y = t), outlier.shape = NA) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  xlab("Cell") + ylab("Pseudotime")

gplt <- plot_grid(mbplt, dplt, labels = c("A", "B"))
cowplot::ggsave(gplt, file = file.path(fig_path, "S2_monocle_comparison.png"), width = 8, height = 3, scale = 1.2)





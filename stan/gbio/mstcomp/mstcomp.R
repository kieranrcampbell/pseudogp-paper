library(monocle)
library(scater)
library(rhdf5)
library(readr)
library(coda)
library(MCMCglmm)

base_dir <- "~/mount/"
setwd(file.path(base_dir, "GP/pseudogp2/stan/gbio/mstcomp"))
source(file.path(base_dir, "GP/pseudogp2/stan/diffexpr/monocle/prep_data.R"))
sce <- load_data(base_dir)

h5file <- file.path(base_dir, "GP/pseudogp2/data/stan_traces_for_gbio.h5")
X <- h5read(h5file, "X")

sce <- plotPCA(sce, return_SCESet = TRUE)
plotReducedDim(sce)

# cds <- toCellDataSet(sce)
# sce2 <- fromCellDataSet(cds, use_exprs_as = "exprs")
# plotReducedDim(sce2)


# comparison of pseudotimes -----------------------------------------------

redDim(sce) <- X
plotReducedDim(sce)
cds <- toCellDataSet(sce)
cds <- orderCells(cds)
cds_le_pseudotime <- cds$Pseudotime

# df <- read_csv(file.path("~/mount", "GP/pseudogp2/R_notebooks/df_fit.csv"))
# to_use <- df$gene_names[df$for_embedding]
# cds <- setOrderingFilter(cds, to_use)
# cds <- reduceDimension(cds, use_irlba = FALSE)
# cds <- orderCells(cds)
# cds_ica_pseudotime <- cds$Pseudotime
# plot_spanning_tree(cds)

tmcmc <- mcmc(h5read(h5file, "pst"))
tmap <- posterior.mode(tmcmc)
intervals <- HPDinterval(tmcmc, prob = 0.95)
dft <- data.frame(monocle = cds_le_pseudotime, tmap = tmap, intervals)

R2 <- signif(cor(dft$monocle, dft$tmap)^2, digits = 2)
R2text <- as.expression(substitute(italic(R)^2 == r, list(r = R2)))

cplt <- ggplot(dft) + geom_point(aes(x = monocle, y = tmap)) + 
  geom_errorbar(aes(x = monocle, ymin = lower, ymax = upper), alpha = 0.6, width = 0.3) +
  xlab("Monocle estimate") + ylab("GPLVM MAP estimate") +
  geom_text(x = 3, y = 0.9, label = as.character(R2text), parse = TRUE)


set.seed(123)
nsamples <- 30
standardize <- function(t) (t - min(t))/(max(t) - min(t))
ncells <- ncol(sce)
pseudosamples <- vector("list", ncells)
gt <- orderCells(cds)$Pseudotime
gt <- standardize(gt)

# did you ever see such awful R in your life
for(i in 1:nsamples) {
  chosen <- sample(ncells, round(0.7*ncells), replace=FALSE)
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

plot(gt, pst_means)

dplt <- ggplot(data.frame(pst_95)) + geom_density(aes(x = pst_95), fill = "darkred") +
  xlab(expression(2 * sigma)) + ylab("Density")

pg <- cowplot::plot_grid(cplt, dplt, labels = c("A", "B"), nrow = 1)
cowplot::ggsave(pg, file = "S2_monocle_comparison.png", width = 8, height = 3, scale = 1.2)


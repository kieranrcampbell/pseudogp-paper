Data

Traces are all stored in pseudogp2/data/monocle_multi_traces.h5 with the following structure:
- multi/ Joint fitting acros all reduced dim representations
multi/X, multi/Y, multi/Z Coordinates of laplacian eigenmaps, PCA & t-SNE respectively

- indv/ Individual fitting for each reduced dim representation separately
indv/g1 Laplacian eigenmaps
indv/g2 PCA
indv/g3 tSNE

Differential expression results (p&q values) are stored in
pseudogp2/data/multi_diffexpr/ - self explanatory
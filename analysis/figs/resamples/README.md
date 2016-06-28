Outline of resampling workflow:
  
`0_create_resamples.R` - this subsamples the cells to create reduced dimensional PCA representations as well as an all-cell PCA representation
`1_fit_gplvm.R` - fits pseudotime curves using GPLVM to the subsampled cells
`2_fit_gplvm_all.R` - fits pseudotime curves using GPLVM to the all-cell representation
`3_de_allcells.R` - robust differential expression across pseudotime traces for all cells
`4_create_robust_sceset.R` - analysis output of (3) and saves the results to SCESet

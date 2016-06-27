Outline of resampling workflow:
  
`0_create_resamples.R` - this subsamples the cells to create reduced dimensional PCA representations as well as an all-cell PCA representation
`1_fit_gplvm.R` - fits pseudotime curves using GPLVM to the subsampled cells
`2_fit_gplvm_all.R` - fits pseudtime curves using GPLVM to the all-cell representation

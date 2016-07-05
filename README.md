# pseudogp-paper

This repo accompanies the paper "Order under uncertainty: on the robustness of differential expression analysis using probabilistic models for pseudotime inference". It requires both the [pseudogp](http://www.github.com/kieranrcampbell/pseudogp) and [switchde](http://www.github.com/kieranrcampbell/switchde) R packages as part of the paper, as well as [scater](http://www.github.com/davismcc/scater).

It is  organised into
* `analysis`: Rmarkdown notebooks for all analysis
* `data`: All MCMC traces, embeddings and gene expression `SCESets` 

In theory navigating to the root directory and running `snakemake` will build all figures for the manuscript (deposited in the `figs` directory). However, this should really be run using a compute cluster (see `run_on_cluster` for an example) as there are many heavy computations.

## Authors

Kieran Campbell & Christopher Yau, July 2016


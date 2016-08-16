"""
Snakefile for reproducible analysis of
Order under uncertainty: on the robustness of differential expression analysis using probabilistic models for pseudotime inference

Kieran Campbell, University of Oxford
July 2015
"""


import glob
import numpy as np

np.random.seed(123)

DE_RUNS = [str(i) for i in range(1, 501)]

BOOTSTRAPS = [i for i in range(1, 501)] # we do 500 bootstraps, but R indexes at 1

RESAMPLES = [str(i) for i in range(1, 101)]

resamples = np.arange(1, 101)
resamples_failed_qc = np.array([2, 3, 7, 11, 20, 26, 41, 52, 55, 59, 82, 85, 86, 87, 94, 96])
resamples_de = np.setdiff1d(resamples, resamples_failed_qc)
resamples_de = np.random.choice(resamples_de, 50, replace = False) # which resampled PCAs do we take for diff expression?
RESAMPLES_DE = [str(i) for i in list(resamples_de)]

TRACE_SAMPLES = [str(i) for i in list(np.random.choice(500, 100, replace = False))]

studies = ["trapnell", "shin", "burns"]

pst_traces = expand("data/{study}_pseudotime_traces.h5", study = studies)
embeddings = expand("data/{study}_embeddings.h5", study = studies)
sces = expand("data/sce_{study}.Rdata", study = studies)

trapnell_de = expand("data/diffexpr/trapnell/de_{tde_run}.csv", tde_run = DE_RUNS)
shin_de = expand("data/diffexpr/shin/de_{sde_run}.csv", sde_run = DE_RUNS)
burns_de = expand("data/diffexpr/burns/de_{bde_run}.csv", bde_run = DE_RUNS)

resample_all_cell_pvals = expand("data/resamples/all_cells_diffexpr/pvals_{posterior_trace}.csv", posterior_trace = TRACE_SAMPLES)


resample_traces = expand("data/resamples/gplvm_fits/fit_{resample}.Rdata", resample = RESAMPLES)
resample_de = expand("data/resamples/diffexpr/pvals_{de_resample}.csv", de_resample = RESAMPLES_DE)

trace_de = expand("data/resamples/trace_diffexpr/pvals_{trace_resample}_{trace}.csv", 
					trace_resample = RESAMPLES_DE, trace = TRACE_SAMPLES)

de_agg = expand("data/diffexpr/agg_pvals/{study}.csv", study = studies)
de_map = expand("data/diffexpr/map/{map_study}.csv", map_study = studies)
fdr_csv = expand("data/diffexpr/{fdr_study}_fdr.csv", fdr_study = studies)

bootstrap_de = expand("data/bootstrap/bootstrapped_de/de_{bs}.csv", bs = BOOTSTRAPS)
gplvm_bootstrap_de = expand("data/bootstrap/gplvm_de/de_{bs}.csv", bs = BOOTSTRAPS)

rule all:
	input:
		pst_traces,
		sces,
		"figs/chains/chains.png",
		"figs/envelope/2_cloud.png",
		"figs/postuncert/pu_density.png",
		"figs/postuncert/3bcd_post_uncert.png",
		"figs/switchres/trapnell_5_switchres.png",
		"figs/switchres/burns_5_switchres.png",
		"figs/switchres/shin_5_switchres.png",
		trapnell_de, shin_de, burns_de,
		resample_all_cell_pvals,
		de_agg, de_map,
		trace_de,
		"figs/diffexpr/go_enriched_categories.png",
		bootstrap_de,
		gplvm_bootstrap_de

"""
		resample_traces,
		resample_de, trace_de,
		de_agg, de_map,
		fdr_csv
"""	

rule trapnell_basic:
	output:
		"data/trapnell_embeddings.h5", # reduced dimension representation
		"data/sce_trapnell.Rdata", # Scater object holding gene expression
		"data/trapnell_pseudotime_traces.h5", # pseudotime traces
		"figs/diagnostic/trapnell.pdf"
	shell:
		"Rscript analysis/basic/trapnell.R"


rule shin_basic:
	output:
		"data/shin_embeddings.h5", # reduced dimension representation
		"data/sce_shin.Rdata", # Scater object holding gene expression
		"data/shin_pseudotime_traces.h5", # pseudotime traces
		"data/waterfall_data.xlsx", # raw
		"figs/diagnostic/shin.pdf"
	shell:
		"Rscript analysis/basic/shin.R"


rule burns_basic:
	output:
		"data/burns_embeddings.h5", # reduced dimension representation
		"data/sce_burns.Rdata", # Scater object holding gene expression
		"data/burns_pseudotime_traces.h5", # pseudotime traces
		"data/burns_raw.txt.gz", # raw TPM data,
		"figs/diagnostic/burns.pdf"
	shell:
		"Rscript analysis/basic/burns.R"

rule fig_chains:
	input:
		"data/trapnell_embeddings.h5"
	output:
		"figs/chains/chains.png",
		"data/chains/fits.Rdata",
	shell:
		"Rscript analysis/figs/chains/chains.R"

rule fig_envelope:
	input:
		pst_traces,
		embeddings,
	output:
		"figs/envelope/2_cloud.png"
	shell:
		"Rscript analysis/figs/envelope/all_uncertainty_paper.R"

rule fig_postuncert:
	input:
		pst_traces
	output:
		"figs/postuncert/pu_density.png",
		"figs/postuncert/3bcd_post_uncert.png"
	shell:
		"Rscript analysis/figs/postuncert/postuncert.R"

rule fig_switchres:
	input:
		sces,
		pst_traces
	output:
		"figs/switchres/trapnell_5_switchres.png",
		"figs/switchres/burns_5_switchres.png",
		"figs/switchres/shin_5_switchres.png"
	shell:
		"Rscript analysis/figs/switchres/switchres.R"


rule trapnell_diffexpr:
	input:
		sce = "data/sce_trapnell.Rdata",
		traces = "data/trapnell_pseudotime_traces.h5",
	output:
		"data/diffexpr/trapnell/de_{tde_run}.csv"
	shell:
		"Rscript analysis/diffexpr/0_diffexpr_analysis.R {input.traces} {input.sce} {wildcards.tde_run} {output}"


rule shin_diffexpr:
	input:
		sce = "data/sce_shin.Rdata",
		traces = "data/shin_pseudotime_traces.h5"
	output:
		"data/diffexpr/shin/de_{sde_run}.csv"
	shell:
		"Rscript analysis/diffexpr/0_diffexpr_analysis.R {input.traces} {input.sce} {wildcards.sde_run} {output}"



rule burns_diffexpr:
	input:
		sce = "data/sce_burns.Rdata",
		traces = "data/burns_pseudotime_traces.h5",
	output:
		"data/diffexpr/burns/de_{bde_run}.csv"
	shell:
		"Rscript analysis/diffexpr/0_diffexpr_analysis.R {input.traces} {input.sce} {wildcards.bde_run} {output}"

rule aggregate_de:
	input:
		trapnell_de, shin_de, burns_de
	output:
		"data/diffexpr/agg_pvals/{study}.csv"
	shell:
		"Rscript analysis/diffexpr/1_aggregate.R {wildcards.study} {output}"

rule map_de:
	input:
		sces, pst_traces
	output:
		"data/diffexpr/map/{map_study}.csv"
	shell:
		"Rscript analysis/diffexpr/2_map_estimates.R {wildcards.map_study} {output}"

rule fdr_calc:
	input:
		de_agg, de_map
	output:
		"data/diffexpr/all_pvals.csv",
		"figs/diffexpr/de_gene_nums.png"
	shell:
		"Rscript analysis/diffexpr/3_fdr_calculation.R"

rule go_analysis:
	input:
		"data/diffexpr/all_pvals.csv"
	output:
		"data/diffexpr/go_no_direction.csv"
	shell:
		"Rscript analysis/diffexpr/4_go_analysis.R"

rule go_figs:
	input:
		"data/diffexpr/go_no_direction.csv"
	output:
		"figs/diffexpr/go_enriched_categories.png",
		"figs/diffexpr/venn_Burns.png",
		"figs/diffexpr/venn_Shin.png",
		"figs/diffexpr/venn_Trapnell.png"
	shell:
		"Rscript analysis/diffexpr/5_go_figs.R"


# Create PCA representations of 80% subsamples & for all cells
rule resample_pca:
	input:
		"data/sce_trapnell.Rdata"
	output:
		"data/resamples/pca_resamples.Rdata", "data/resamples/pca_all.Rdata"
	shell:
		"Rscript analysis/figs/resamples/0_create_resamples.R"

# Fit GPLVMs for each PCA resample
rule resample_gplvm:
	input:
		"data/resamples/pca_resamples.Rdata"
	output:
		"data/resamples/gplvm_fits/fit_{resample}.Rdata"
	shell:
		"Rscript analysis/figs/resamples/1_fit_gplvm.R {wildcards.resample}"

# Fit GPLVM for all-cell PCA
rule resample_gplvm_all:
	input:
		"data/resamples/pca_all.Rdata"
	output:
		"data/resamples/gplvm_fit_all.Rdata"
	shell:
		"Rscript analysis/figs/resamples/2_fit_gplvm_all.R"


# Perform DE across traces for the all-cell PCA
rule resample_all_diffexpr:
	input:
		"data/resamples/gplvm_fit_all.Rdata", "data/sce_trapnell.Rdata",
		"data/resamples/pca_all.Rdata"
	output:
		"data/resamples/all_cells_diffexpr/pvals_{posterior_trace}.csv"
	shell:
		"Rscript analysis/figs/resamples/3_de_allcells.R {wildcards.posterior_trace}"

# Use only robustly DE genes at all-cell estimate for subsample differential expression
rule resample_create_robust_sce:
	input:
		resample_all_cell_pvals
	output:
		"data/resamples/sce_trapnell_robust.Rdata"
	shell:
		"Rscript analysis/figs/resamples/4_create_robust_sceset.R"

# Differential expression across all subsamples
rule trace_diffexpr:
	input:
		"data/resamples/sce_trapnell_robust.Rdata",
		"data/resamples/gplvm_fits/fit_{trace_resample}.Rdata",
		"data/resamples/pca_resamples.Rdata"
	output:
		"data/resamples/trace_diffexpr/pvals_{trace_resample}_{trace}.csv"
	shell:
		"Rscript analysis/figs/resamples/5_trace_de.R {wildcards.trace_resample} {wildcards.trace}"

## Bootstrapping -------

rule create_bootstraps:
	input:
		"data/sce_trapnell.Rdata"
	output:
		"data/bootstrap/bootstrap_pseudotimes.csv",
		"data/bootstrap/which_cells.csv"
	shell:
		"Rscript analysis/figs/boostrap/0_create_bootstraps.R"

rule create_gplvm_pst:
	input:
		"data/sce_trapnell.Rdata"
	output:
		"data/bootstrap/gplvm_pseudotimes.csv"
	shell:
		"Rscript analysis/figs/boostrap/1_fit_gplvm.R"

rule bootstrap_de:
	input:
		"data/sce_trapnell.Rdata",
		"data/bootstrap/which_cells.csv",
		"data/bootstrap/bootstrap_pseudotimes.csv"
	output:
		"data/bootstrap/bootstrapped_de/de_{bs}.csv"
	shell:
		"Rscript analysis/figs/boostrap/2_bootstrap_de.R {bs}"

rule bootstrap_gplvm_de:
	input:
		"data/sce_trapnell.Rdata",
		"data/bootstrap/gplvm_pseudotimes.csv"
	output:
		"data/bootstrap/gplvm_de/de_{bs}.csv"
	shell:
		"Rscript analysis/figs/boostrap/3_gplmv_de.R {bs}"



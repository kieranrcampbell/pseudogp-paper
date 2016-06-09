import glob
import numpy as np

np.random.seed(123)

#DE_RUNS = ["1x50", "51x100", "101x150", "151x200", "201x250", "251x300",
#			"301x350", "351x400", "401x450", "451x500"]
DE_RUNS = ["1x4", "5x9"]

RESAMPLES = [str(i) for i in range(1, 101)]

resamples = np.arange(1, 101)
resamples_failed_qc = np.array([2, 3, 7, 11, 20, 26, 41, 52, 55, 59, 82, 85, 86, 87, 94, 96])
resamples_de = np.setdiff1d(resamples, resamples_failed_qc)
resamples_de = np.random.choice(resamples_de, 50)
RESAMPLES_DE = [str(i) for i in list(resamples_de)]

studies = ["trapnell", "shin", "burns"]

pst_traces = expand("data/{study}_pseudotime_traces.h5", study = studies)
embeddings = expand("data/{study}_embeddings.h5", study = studies)
sces = expand("data/sce_{study}.Rdata", study = studies)

trapnell_de = expand("data/diffexpr/trapnell/de_{tde_run}.csv", tde_run = DE_RUNS)
shin_de = expand("data/diffexpr/shin/de_{sde_run}.csv", sde_run = DE_RUNS)
burns_de = expand("data/diffexpr/burns/de_{bde_run}.csv", bde_run = DE_RUNS)

resample_traces = expand("data/resamples/gplvm_fits/fit_{resample}.Rdata", resample = RESAMPLES)
resample_de = expand("data/resamples/diffexpr/pvals_{de_resample}.csv", de_resample = RESAMPLES_DE)

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
		#"figs/fdr.png",
		trapnell_de, shin_de, burns_de,
		resample_traces,
		resample_de
	

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
		"Rscript analysis/diffexpr/0_diffexpr_analysis.R {input.traces} {input.sce} {output}"

rule trapnell_fdr:
	input:
		traces = "data/trapnell_pseudotime_traces.h5",
		de = trapnell_de,
		sce = "data/sce_trapnell.Rdata"
	output:
		plot_pdf = "figs/diffexpr/trapnell_plots.pdf",
		fdr_csv = "data/diffexpr/trapnell_fdr.csv"
	shell:
		"Rscript analysis/diffexpr/1_fdr_calculation.R {input.traces} {input.de} {input.sce} {output.plot_pdf} {output.fdr_csv}"

rule shin_diffexpr:
	input:
		sce = "data/sce_shin.Rdata",
		traces = "data/shin_pseudotime_traces.h5"
	output:
		"data/diffexpr/shin/de_{sde_run}.csv"
	shell:
		"Rscript analysis/diffexpr/0_diffexpr_analysis.R {input.traces} {input.sce} {output}"


rule shin_fdr:
	input:
		traces = "data/shin_pseudotime_traces.h5",
		de = shin_de,
		sce = "data/sce_shin.Rdata"
	output:
		plot_pdf = "figs/diffexpr/shin_plots.pdf",
		fdr_csv = "data/diffexpr/shin_fdr.csv"
	shell:
		"Rscript analysis/diffexpr/1_fdr_calculation.R {input.traces} {input.de} {input.sce} {output.plot_pdf} {output.fdr_csv}"


rule burns_diffexpr:
	input:
		sce = "data/sce_burns.Rdata",
		traces = "data/burns_pseudotime_traces.h5",
	output:
		"data/diffexpr/burns/de_{bde_run}.csv"
	shell:
		"Rscript analysis/diffexpr/0_diffexpr_analysis.R {input.traces} {input.sce} {output}"

rule burns_fdr:
	input:
		traces = "data/burns_pseudotime_traces.h5",
		de = burns_de,
		sce = "data/sce_burns.Rdata"
	output:
		plot_pdf = "figs/diffexpr/burns_plots.pdf",
		fdr_csv = "data/diffexpr/burns_fdr.csv"
	shell:
		"Rscript analysis/diffexpr/1_fdr_calculation.R {input.traces} {input.de} {input.sce} {output.plot_pdf} {output.fdr_csv}"

rule make_fdr_plots:
	input:
		"data/diffexpr/burns_fdr.csv",
		"data/diffexpr/trapnell_fdr.csv",
		"data/diffexpr/shin_fdr.csv",
		trapnell_de, burns_de, shin_de,
		"data/sce_trapnell.Rdata"
	output:
		"figs/fdr.png"
	shell:
		"Rscript analysis/diffexpr/fdr_plots.R"

rule resample_pca:
	input:
		"data/sce_trapnell.Rdata"
	output:
		"data/resamples/pca_resamples.Rdata"
	shell:
		"Rscript analysis/figs/resamples/0_create_resamples.R"

rule resample_gplvm:
	input:
		"data/resamples/pca_resamples.Rdata"
	output:
		"data/resamples/gplvm_fits/fit_{resample}.Rdata"
	shell:
		"Rscript analysis/figs/resamples/1_fit_gplvm.R {wildcards.resample}"

rule resample_choose_genes:
	input:
		"data/sce_trapnell.Rdata"
	output:
		"data/resamples/sce_trapnell_resamples.Rdata"
	shell:
		"Rscript analysis/figs/resamples/2_choose_genes.R"

rule resample_diffexpr:
	input:
		"data/resamples/sce_trapnell_resamples.Rdata",
		"data/resamples/gplvm_fits/fit_{de_resample}.Rdata"
	output:
		"data/resamples/diffexpr/pvals_{de_resample}.csv"
	shell:
		"Rscript analysis/figs/resamples/3_differential_expression.R {wildcards.de_resample}"


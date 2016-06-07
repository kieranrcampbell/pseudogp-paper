import glob

de_runs = ["1x50", "51x100", "101x150", "151x200", "201x250", "251x300",
			"301x350", "351x400", "401x450", "451x500"]

studies = ["trapnell", "shin", "burns"]

pst_traces = expand("data/{study}_pseudotime_traces.h5", study = studies)
embeddings = expand("data/{study}_embeddings.h5", study = studies)
sces = expand("data/sce_{study}.Rdata", study = studies)


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
		"figs/fdr.png"

	

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

rule tmp_confuse_snakemake:
	output:
		expand("data/diffexpr/trapnell_{de_run}", de_run = de_runs),
		expand("data/diffexpr/shin_{de_run}", de_run = de_runs),
		expand("data/diffexpr/burns_{de_run}", de_run = de_runs)
	run:
		for de_run in de_runs:
		    shell("touch data/diffexpr/trapnell_{de_run} data/diffexpr/shin_{de_run} data/diffexpr/burns_{de_run}")

rule trapnell_diffexpr:
	input:
		sce = "data/sce_trapnell.Rdata",
		traces = "data/trapnell_pseudotime_traces.h5",
		tmp = expand("data/diffexpr/trapnell_{de_run}", de_run = de_runs)
	output:
		traces = "data/diffexpr/trapnell_de_traces.h5",
	shell:
		"Rscript analysis/diffexpr/0_diffexpr_analysis.R {input.traces} {input.sce} {output.traces} {input.tmp}"

rule trapnell_fdr:
	input:
		traces = "data/trapnell_pseudotime_traces.h5",
		de = "data/diffexpr/trapnell_de_traces.h5",
		sce = "data/sce_trapnell.Rdata"
	output:
		plot_pdf = "figs/diffexpr/trapnell_plots.pdf",
		fdr_csv = "data/diffexpr/trapnell_fdr.csv"
	shell:
		"Rscript analysis/diffexpr/1_fdr_calculation.R {input.traces} {input.de} {input.sce} {output.plot_pdf} {output.fdr_csv}"

rule shin_diffexpr:
	input:
		sce = "data/sce_shin.Rdata",
		traces = "data/shin_pseudotime_traces.h5",
		tmp = expand("data/diffexpr/shin_{de_run}", de_run = de_runs)
	output:
		traces = "data/diffexpr/shin_de_traces.h5",
	shell:
		"Rscript analysis/diffexpr/0_diffexpr_analysis.R {input.traces} {input.sce} {output.traces} {input.tmp}"


rule shin_fdr:
	input:
		traces = "data/shin_pseudotime_traces.h5",
		de = "data/diffexpr/shin_de_traces.h5",
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
		tmp = expand("data/diffexpr/burns_{de_run}", de_run = de_runs)
	output:
		traces = "data/diffexpr/burns_de_traces.h5",
	shell:
		"Rscript analysis/diffexpr/0_diffexpr_analysis.R {input.traces} {input.sce} {output.traces} {input.tmp}"

rule burns_fdr:
	input:
		traces = "data/burns_pseudotime_traces.h5",
		de = "data/diffexpr/burns_de_traces.h5",
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
		"data/diffexpr/trapnell_de_traces.h5",
		"data/trapnell_pseudotime_traces.h5",
		"data/sce_trapnell.Rdata"
	output:
		"figs/fdr.png"
	shell:
		"Rscript analysis/diffexpr/fdr_plots.R"





library(ggplot2)
library(tidyr)
library(magrittr)
library(goseq)
library(readr)
library(dplyr)
library(rhdf5)
library(MCMCglmm)

rename <- dplyr::rename
select <- dplyr::select
mutate <- dplyr::mutate
filter <- dplyr::filter



alpha <- 0.05


do_go_enrichment <- function(t_df, sce, pseudotime, gene_directions, genome, id) {
  
  t_df <- mutate(t_df, signif_map = q_val < alpha, 
                   signif_gplvm = prop_sig > (1 - alpha),
                 direction = gene_directions)
  t_df <- mutate(t_df, up_signif_map = signif_map & direction == 1,
                 up_signif_gplvm = signif_gplvm & direction == 1,
                 down_signif_map = signif_map & direction == -1,
                 down_signif_gplvm = signif_gplvm & direction == -1)
  
  genes <- data.frame(up_map = as.integer(t_df$up_signif_map),
                up_gplvm = as.integer(t_df$up_signif_gplvm),
                down_map = as.integer(t_df$down_signif_map),
                down_gplvm = as.integer(t_df$down_signif_gplvm))
  rownames(genes) <- t_df$gene
  
  do_goseq <- function(genes, genome = "hg19", id = "ensGene", test.cats = "GO:BP") {
    pwf <- nullp(genes, genome, id)
    go <- goseq(pwf, genome, id, test.cats = test.cats)
    go <- mutate(go, qval = p.adjust(over_represented_pvalue, method = "BH"))
    return(go)
  }
  
  go_results <- apply(genes, 2, do_goseq, genome, id)
  
  up <- inner_join(go_results$up_map, go_results$up_gplvm, 
                  by = c("category", "term", "ontology"))  %>%
    dplyr::select(-starts_with("over_represented"), -starts_with("under_represented"), -starts_with("num")) %>%
    rename(qval_map = qval.x, qval_gplvm = qval.y) %>%
    mutate(signif_map = qval_map < alpha, signif_gplvm = qval_gplvm < alpha) %>% tbl_df()
  
  down <- inner_join(go_results$down_map, go_results$down_gplvm, 
                   by = c("category", "term", "ontology"))  %>%
    dplyr::select(-starts_with("over_represented"), -starts_with("under_represented"), -starts_with("num")) %>%
    rename(qval_map = qval.x, qval_gplvm = qval.y) %>%
    mutate(signif_map = qval_map < alpha, signif_gplvm = qval_gplvm < alpha) %>% tbl_df()
  
  up <- mutate(up, direction = "up")
  down <- mutate(down, direction = "down")
  
  results <- rbind(up, down)
  return(results)
}

sig_df <- read_csv("data/diffexpr/all_pvals.csv") # data frame for all results

# Trapnell ----------------------------------------------------------------

parse_trapnell_gene_names <- function(g) {
  sapply(strsplit(g, ".", fixed = TRUE), `[[`, 1)
}

t_df <- filter(sig_df, study == "trapnell")


load("data/sce_trapnell.Rdata")
pst <- h5read("data/trapnell_pseudotime_traces.h5", "pst")
pseudotime <- posterior.mode(mcmc(pst))

gex_vector <- exprs(sce[t_df$gene, ])
trapnell_gene_directions <- apply(gex_vector, 1, function(x) {
  sign(cor(x, pseudotime))
})

t_df <- mutate(t_df, gene = parse_trapnell_gene_names(gene))
go_trapnell <- do_go_enrichment(t_df, sce, pseudotime, trapnell_gene_directions, "hg19", "ensGene")



# Shin --------------------------------------------------------------------
s_df <- filter(sig_df, study == "shin")
load("data/sce_shin.Rdata")
pst <- h5read("data/shin_pseudotime_traces.h5", "pst")
pseudotime <- posterior.mode(mcmc(pst))

gex_vector <- exprs(sce[s_df$gene, ])
shin_gene_directions <- apply(gex_vector, 1, function(x) {
  sign(cor(x, pseudotime))
})
go_shin <- do_go_enrichment(s_df, sce, pseudotime, shin_gene_directions, "mm10", "geneSymbol")


# Burns -------------------------------------------------------------------

b_df <- filter(sig_df, study == "burns")
load("data/sce_burns.Rdata")
pst <- h5read("data/burns_pseudotime_traces.h5", "pst")
pseudotime <- posterior.mode(mcmc(pst))

gex_vector <- exprs(sce[b_df$gene, ])
burns_gene_directions <- apply(gex_vector, 1, function(x) {
  sign(cor(x, pseudotime))
})

go_burns <- do_go_enrichment(b_df, sce, pseudotime, burns_gene_directions, "mm10", "geneSymbol")



# Analyse together --------------------------------------------------------

go_trapnell <- mutate(go_trapnell, study = "trapnell")
go_shin <- mutate(go_shin, study = "shin")
go_burns <- mutate(go_burns, study = "burns")

go <- rbind(go_trapnell, go_shin, go_burns)

write_csv(go, "data/diffexpr/go.csv")



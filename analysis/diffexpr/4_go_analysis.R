
library(ggplot2)
library(tidyr)
library(magrittr)
library(goseq)
library(readr)
library(dplyr)

rename <- dplyr::rename
select <- dplyr::select
mutate <- dplyr::mutate

sig_df <- read_csv("data/diffexpr/all_pvals.csv")


do_go_enrichment <- function(t_df, genome, id) {
  t_df <- mutate(t_df, signif_map = q_val < alpha, 
                 signif_gplvm = prop_sig > (1 - alpha))
  t_df <- mutate(t_df, map_only = signif_map & !signif_gplvm)

  genes <- data.frame(map = as.integer(t_df$signif_map),
                      gplvm = as.integer(t_df$signif_gplvm),
                      map_only = as.integer(t_df$map_only))
  rownames(genes) <- t_df$gene
  
  do_goseq <- function(genes, genome = "hg19", id = "ensGene", test.cats = "GO:BP") {
    pwf <- nullp(genes, genome, id)
    go <- goseq(pwf, genome, id, test.cats = test.cats)
    go <- mutate(go, qval = p.adjust(over_represented_pvalue, method = "BH"))
    return(go)
  }
  
  go_results <- apply(genes, 2, do_goseq, genome, id)
  gres <- Map(function(df, name) {
    mutate(df, etype = name)
  }, go_results, names(genes))
  gres <- do.call(rbind, gres) %>% tbl_df()
  
  return(gres)
}



# Trapnell ----------------------------------------------------------------
t_df <- filter(sig_df, study == "trapnell")

alpha <- 0.05

parse_trapnell_gene_names <- function(g) {
  sapply(strsplit(g, ".", fixed = TRUE), `[[`, 1)
}
t_df <- mutate(t_df, gene = parse_trapnell_gene_names(gene))

go_trapnell <- do_go_enrichment(t_df, "hg19", "ensGene")

# Shin --------------------------------------------------------------------
s_df <- filter(sig_df, study == "shin")
go_shin <- do_go_enrichment(s_df, "mm10", "geneSymbol")


# Burns -------------------------------------------------------------------

b_df <- filter(sig_df, study == "burns")
go_burns <- do_go_enrichment(b_df, "mm10", "geneSymbol")


# Save --------------------------------------------------------------------

go_trapnell <- mutate(go_trapnell, study = "trapnell")
go_shin <- mutate(go_shin, study = "shin")
go_burns <- mutate(go_burns, study = "burns")

go <- rbind(go_trapnell, go_shin, go_burns)

write_csv(go, "data/diffexpr/go_no_direction.csv")
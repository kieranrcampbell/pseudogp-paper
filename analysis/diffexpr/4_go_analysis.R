library(dplyr)
library(ggplot2)
library(tidyr)
library(magrittr)
library(goseq)

sig_df <- read_csv("data/diffexpr/all_pvals.csv")


# Trapnell ----------------------------------------------------------------
t_df <- filter(sig_df, study == "trapnell")

alpha <- 0.05

parse_trapnell_gene_names <- function(g) {
  sapply(strsplit(g, ".", fixed = TRUE), `[[`, 1)
}
t_df <- mutate(t_df, gene = parse_trapnell_gene_names(gene))

t_df <- mutate(t_df, signif_map = q_val < alpha, 
                 signif_gplvm = prop_sig > (1 - alpha))

genes_map <- as.integer(t_df$signif_map)
genes_gplvm <- as.integer(t_df$signif_gplvm)
genes_random <- rep(0, length(genes_map))
genes_random[sample(which(genes_map == 1), sum(genes_gplvm == 1))] <- 1
names(genes_map) <- names(genes_gplvm) <- names(genes_random) <- t_df$gene

do_goseq <- function(genes, genome = "hg19", id = "ensGene", test.cats = "GO:BP") {
  pwf <- nullp(genes, genome, id)
  go <- goseq(pwf, genome, id, test.cats = test.cats)
  go <- mutate(go, qval = p.adjust(over_represented_pvalue, method = "BH"))
  return(go)
}

go_map <- do_goseq(genes_map)
go_gplvm <- do_goseq(genes_gplvm)
go_random <- do_goseq(genes_random)

## number significant in each case
sum(go_map$qval < alpha)
sum(go_gplvm$qval < alpha)
sum(go_random$qval < alpha)

a <- tbl_df(inner_join(go_map, go_gplvm, 
                       by = c("category", "term", "ontology"))) %>%
  select(-starts_with("over_represented"), 
         -starts_with("under_represented"), 
         -starts_with("num")) %>%
  rename(qval_map = qval.x, qval_gplvm = qval.y) %>%
  mutate(signif_map = qval_map < alpha, signif_gplvm = qval_gplvm < alpha)

ggplot(a, aes(x = qval_map, y = qval_gplvm, colour = signif_map)) + 
  geom_point(alpha = 0.5) + scale_color_brewer(palette = "Set1")

a_map <- filter(a, qval_gplvm > 0.05, qval_map < alpha) %>% arrange(qval_map) #%>% View()
a_gp <- filter(a, qval_map > 0.05, qval_gplvm < alpha) %>% arrange(qval_gplvm) 

a_gp$term[grep("mito", a_gp$term)]
a_map$term[grep("mito", a_map$term)]

a_gp$term[grep("muscle", a_gp$term)]
a_map$term[grep("muscle", a_map$term)]

a_gp$term[grep("cell cycle", a_gp$term)]
a_map$term[grep("cell cycle", a_map$term)]

a %>% filter(qval_gplvm < 0.01) %>% select(category, qval_gplvm) %>%
  write_delim("~/Desktop/gp.txt")
a %>% filter(qval_map < 0.01) %>% select(category, qval_map) %>%
  write_delim("~/Desktop/map.txt")

muscle <- a[grep("mito", a$term),] %>%
  gather(measure, qval, qval_map, qval_gplvm) %>%
  filter(qval < 0.01)
ggplot(muscle, aes(x = measure, y = log10(qval))) + geom_violin()


# Shin --------------------------------------------------------------------
s_df <- filter(sig_df, study == "shin")

s_df <- mutate(s_df, signif_map = q_val < alpha, 
               signif_gplvm = prop_sig > (1 - alpha))

shin_genes_map <- as.integer(s_df$signif_map)
shin_genes_gplvm <- as.integer(s_df$signif_gplvm)
names(shin_genes_map) <- names(shin_genes_gplvm) <- s_df$gene

shin_pwf_map <- nullp(shin_genes_map, "mm10", "geneSymbol")
shin_go_map <- goseq(shin_pwf_map, "mm10", "geneSymbol", test.cats = "GO:BP")
shin_go_map <- mutate(go_map, qval = p.adjust(over_represented_pvalue, method = "BH"))

shin_pwf_gplvm <- nullp(shin_genes_gplvm, "mm10", "geneSymbol")
shin_go_gplvm <- goseq(shin_pwf_gplvm, "mm10", "geneSymbol", test.cats = "GO:BP")
shin_go_gplvm <- mutate(shin_go_gplvm, qval = p.adjust(over_represented_pvalue, method = "BH"))


ashin <- tbl_df(inner_join(shin_go_map, shin_go_gplvm, by = c("category", "term", "ontology"))) %>%
  select(-starts_with("over_represented"), -starts_with("under_represented"), -starts_with("num")) %>%
  rename(qval_map = qval.x, qval_gplvm = qval.y) %>%
  mutate(signif_map = qval_map < alpha, signif_gplvm = qval_gplvm < alpha)

ggplot(ashin, aes(x = qval_map, y = qval_gplvm, colour = signif_map)) + 
  geom_point(alpha = 0.5)

#a_map <- filter(a, qval_gplvm > 0.05, qval_map < alpha) %>% arrange(qval_map) #%>% View()
#a_gp <- filter(a, qval_map > 0.05, qval_gplvm < alpha) %>% arrange(qval_gplvm) 


# Burns -------------------------------------------------------------------

b_df <- filter(sig_df, study == "burns")
b_df <- mutate(b_df, signif_map = q_val < alpha, 
               signif_gplvm = prop_sig > (1 - alpha))

burns_genes_map <- as.integer(b_df$signif_map)
burns_genes_gplvm <- as.integer(b_df$signif_gplvm)
names(burns_genes_map) <- names(burns_genes_gplvm) <- b_df$gene

burns_pwf_map <- nullp(burns_genes_map, "mm10", "geneSymbol")
burns_go_map <- goseq(burns_pwf_map, "mm10", "geneSymbol", test.cats = "GO:BP")
burns_go_map <- mutate(burns_go_map, qval = p.adjust(over_represented_pvalue, method = "BH"))

burns_pwf_gplvm <- nullp(burns_genes_gplvm, "mm10", "geneSymbol")
burns_go_gplvm <- goseq(burns_pwf_gplvm, "mm10", "geneSymbol", test.cats = "GO:BP")
burns_go_gplvm <- mutate(burns_go_gplvm, qval = p.adjust(over_represented_pvalue, method = "BH"))

aburns <- tbl_df(inner_join(burns_go_map, burns_go_gplvm, by = c("category", "term", "ontology"))) %>%
  select(-starts_with("over_represented"), -starts_with("under_represented"), -starts_with("num")) %>%
  rename(qval_map = qval.x, qval_gplvm = qval.y) %>%
  mutate(signif_map = qval_map < alpha, signif_gplvm = qval_gplvm < alpha)

ggplot(aburns, aes(x = qval_map, y = qval_gplvm, colour = signif_map)) + 
  geom_point(alpha = 0.5)

# a_map <- filter(a, qval_gplvm > 0.05, qval_map < alpha) %>% arrange(qval_map) #%>% View()
# a_gp <- filter(a, qval_map > 0.05, qval_gplvm < alpha) %>% arrange(qval_gplvm) 

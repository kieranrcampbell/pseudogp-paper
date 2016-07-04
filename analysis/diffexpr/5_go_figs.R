
library(ggplot2)
library(dplyr)
library(cowplot)
library(Hmisc)
library(RColorBrewer)
library(VennDiagram)
library(readr)
library(tidyr)

alpha <- 0.05

# Figures -----------------------------------------------------------------


go <- read_csv("data/diffexpr/go_no_direction.csv")
# go <- mutate(go, type = paste(as.character(signif_map), as.character(signif_gplvm), sep = "_"))
# go <- mutate(go, type = plyr::mapvalues(type, from = c("TRUE_TRUE", "TRUE_FALSE", "FALSE_TRUE", "FALSE_FALSE"),
#                                         to = c("Common", "MAP only", "GPLVM only", "Neither")))

go <- mutate(go, study = capitalize(study))
go$etype <- as.factor(go$etype)

go_sum <- filter(go, qval < alpha) %>%
  group_by(study, etype) %>% 
  summarise(count = n()) %>% 
  complete(study, etype) %>% 
  mutate(count, count = replace(count, is.na(count), 0))
  
  ggplot(go_sum, aes(x = study, fill = etype, y = count)) + 
    geom_bar(stat = "identity", position = "dodge", color = 'black') +
    scale_fill_brewer(palette = "Set1",
                      breaks = c("gplvm", "map", "map_only"),
                      labels = c("Robust",
                                 "All",
                                 "Unstable")) +
    cowplot::theme_cowplot() +
    theme(legend.position = "right",
          legend.title = element_blank(),
          legend.text = element_text(size = 16),
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 14)) +
    xlab("Study") + ylab("Number of enriched GO categories") 

ggsave("figs/diffexpr/go_enriched_categories.png", width = 6, height = 4)


cols <- brewer.pal(3, "Set1")

studies <- unique(go$study)

go_lists <- lapply(studies, function(Study) {
  t <- filter(go, study == Study, qval < alpha)
  labels <- c("All",
             "Robust",
             "Unstable")
  map <- filter(t, etype == "map")$category
  gplvm <- filter(t, etype == "gplvm")$category
  map_only <- filter(t, etype == "map_only")$category
  
  ulist <- list(map = map, gplvm = gplvm, map_only = map_only)
  names(ulist) <- labels
  ulist
})

names(go_lists) <- studies

# Need to do each separately

for(study in studies) {
  venn.diagram(go_lists[[study]],
               filename = paste0("figs/diffexpr/venn_", study, ".png"),
               main = study,
               alpha = 0.3, 
               fill = cols, 
               cex = 1.8,
               cat.cex = 1.8,
               cat.default.pos = "outer", margin = 0.1,
               cat.dist = c(0.1, 0.1, 0.05),
               main.cex = 2.1,
               cat.fontface = "bold",
               fontface = "bold",
               main.fontface = "bold",
               scaled = FALSE,
               euler.d = FALSE,
               main.pos = c(0.5, 0.95))
}





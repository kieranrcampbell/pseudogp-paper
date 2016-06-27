
library(ggplot2)
library(dplyr)
library(cowplot)
library(Hmisc)


alpha <- 0.05

# Figures -----------------------------------------------------------------


go <- read_csv("data/diffexpr/go_no_direction.csv")
# go <- mutate(go, type = paste(as.character(signif_map), as.character(signif_gplvm), sep = "_"))
# go <- mutate(go, type = plyr::mapvalues(type, from = c("TRUE_TRUE", "TRUE_FALSE", "FALSE_TRUE", "FALSE_FALSE"),
#                                         to = c("Common", "MAP only", "GPLVM only", "Neither")))

go <- mutate(go, study = capitalize(study))
go$etype <- as.factor(go$etype)

filter(go, qval < alpha) %>%
  ggplot(aes(x = study, fill = etype)) + geom_bar(position = "dodge") +
  scale_fill_brewer(palette = "Set1",
                    name = "Using genes\ndifferentially expressed...",
                    breaks = c("gplvm", "map", "map_only"),
                    labels = c("Across MCMC traces",
                               "Across MCMC traces \nand MAP estimate",
                               "At MAP estimate only")) +
  theme(legend.position = "right",
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 9),
        axis.title.y = element_blank()) +
  xlab("Study") + ggtitle("Number of enriched GO categories")

ggsave("figs/diffexpr/go_enriched_categories.png", width = 6, height = 4)

## upsetr plot

for(Study in unique(go$study)) {
  t <- filter(go, study == Study, qval < alpha)
  labels = c("Across MCMC traces",
             "Across MCMC traces \nand MAP estimate",
             "At MAP estimate only")
  map <- filter(t, etype == "map")$category
  gplvm <- filter(t, etype == "gplvm")$category
  map_only <- filter(t, etype == "map_only")$category
  
  ulist <- list(map = map, gplvm = gplvm, map_only = map_only)
  names(ulist) <- c("Across MCMC traces \nand MAP estimate",
                    "Across MCMC traces",
                    "At MAP estimate only")
  pdf(paste0("figs/diffexpr/go_", Study, ".pdf"), width = 8, height = 6)
  upset(fromList(ulist), order.by = "freq")
        # name.size = 20, point.size = 4, line.size = 1.5) 
  dev.off()
}

# 
# filter(go, type != "Neither") %>%
#   ggplot(aes(x = type, fill = type)) + geom_bar() +
#   facet_wrap(~ study) +
#   scale_fill_brewer(palette = "Set1") +
#   theme(legend.position = "none",
#         axis.text.x = element_text(angle = 45, vjust = 0.8, size = 10, hjust = 0.8),
#         axis.title.x = element_blank()) +
#   ylab("Number of enriched GO terms") 
# 
# ggsave("figs/diffexpr/go_no_direction.png", width=6, height=4)
# 
# 
# # UP-down analysis --------------------------------------------------------
# 
# go <- read_csv("data/diffexpr/go.csv")
# go <- mutate(go, type = paste(as.character(signif_map), as.character(signif_gplvm), sep = "_"))
# go <- mutate(go, type = plyr::mapvalues(type, from = c("TRUE_TRUE", "TRUE_FALSE", "FALSE_TRUE", "FALSE_FALSE"),
#                                       to = c("Common", "MAP only", "GPLVM only", "Neither")))
# 
# go <- mutate(go, study = capitalize(study), direction = capitalize(direction))
# 
# filter(go, type != "Neither") %>%
#   ggplot(aes(x = direction, fill = type)) + geom_bar(position = "dodge") + facet_wrap(~ study) +
#   xlab("Pseudotime regulation") + ylab("Number of enriched GO terms") +
#   scale_fill_brewer(palette = "Set1") +
#   theme(legend.title = element_blank())
# 
# ggsave("figs/diffexpr/go_direction.png", width=6, height=4)
# 

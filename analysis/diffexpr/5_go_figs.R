
library(ggplot2)
library(dplyr)
library(cowplot)
library(Hmisc)



# Figures -----------------------------------------------------------------


go <- read_csv("data/diffexpr/go_no_direction.csv")
go <- mutate(go, type = paste(as.character(signif_map), as.character(signif_gplvm), sep = "_"))
go <- mutate(go, type = plyr::mapvalues(type, from = c("TRUE_TRUE", "TRUE_FALSE", "FALSE_TRUE", "FALSE_FALSE"),
                                        to = c("Common", "MAP only", "GPLVM only", "Neither")))

go <- mutate(go, study = capitalize(study))

filter(go, type != "Neither") %>%
  ggplot(aes(x = type, fill = type)) + geom_bar() +
  facet_wrap(~ study) +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 0.8, size = 10, hjust = 0.8),
        axis.title.x = element_blank()) +
  ylab("Number of enriched GO terms") 

ggsave("figs/diffexpr/go_no_direction.png", width=6, height=4)


# UP-down analysis --------------------------------------------------------

go <- read_csv("data/diffexpr/go.csv")
go <- mutate(go, type = paste(as.character(signif_map), as.character(signif_gplvm), sep = "_"))
go <- mutate(go, type = plyr::mapvalues(type, from = c("TRUE_TRUE", "TRUE_FALSE", "FALSE_TRUE", "FALSE_FALSE"),
                                      to = c("Common", "MAP only", "GPLVM only", "Neither")))

go <- mutate(go, study = capitalize(study), direction = capitalize(direction))

filter(go, type != "Neither") %>%
  ggplot(aes(x = direction, fill = type)) + geom_bar(position = "dodge") + facet_wrap(~ study) +
  xlab("Pseudotime regulation") + ylab("Number of enriched GO terms") +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.title = element_blank())

ggsave("figs/diffexpr/go_direction.png", width=6, height=4)


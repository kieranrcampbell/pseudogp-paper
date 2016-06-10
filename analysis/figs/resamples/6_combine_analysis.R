library(readr)
library(dplyr)
library(matrixStats)
library(reshape2)
library(magrittr)
library(tidyr)

theme_set(theme_bw())

d <- "data/resamples/trace_diffexpr"

files <- dir(d)

resample_qvals <- lapply(files, function(f) {
  split <- strsplit(f, "_")[[1]]
  rep <- as.numeric(split[2])
  trace <- as.numeric(gsub(".csv", "", split[3], fixed = TRUE))
  csv <- read_csv(file.path(d, f))
  csv %<>% mutate(trace = trace, rep = rep)
  return(csv)
})

rq <- do.call(rbind, resample_qvals)

r <- filter(rq, gene == "ENSG00000120798.11")

ggplot(r, aes(x = trace, y = rep, colour = log(q_val))) + geom_point()

by_trace <- group_by(rq, trace) %>%
  summarise(Variance = var((q_val)))

by_rep <- group_by(rq, rep) %>%
  summarise(Variance = var((q_val)))

var_df <- data_frame(DimRed = by_trace$Variance, 
                     GPVLM = by_rep$Variance) %>%
  gather(Source, Variance)
  
ggplot(var_df, aes(x = Source, y = Variance)) + geom_boxplot()

across_reps <- group_by(rq, gene, rep) %>%
  summarise(mean_qval = mean(q_val)) %>%
  summarise(var_qval = var(mean_qval))

across_gplvm <- group_by(rq, gene, rep) %>%
  summarise(var_qval = var(q_val)) %>%
  summarise(mean_var_qval = mean(var_qval))

vdf <- data_frame(gplvm = across_gplvm$mean_var_qval, 
                  dimrep = across_reps$var_qval) %>%
  gather(source, variance)

mean_v <- vdf %>% group_by(source) %>% summarise(m = mean(variance))

ggplot(vdf, aes(x = source, y = variance)) + geom_boxplot()

s <- replicate(100, {
  x <- matrix(rnorm(1000), nrow = 100, ncol = 10)
  c(mean(rowVars(x)), var(rowMeans(x)))
})

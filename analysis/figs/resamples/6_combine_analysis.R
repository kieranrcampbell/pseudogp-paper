library(scater)
library(embeddr)
library(MCMCglmm)
library(coda)
library(rstan)
library(readr)
library(dplyr)
library(matrixStats)
library(reshape2)

theme_set(theme_bw())

resample_dir <- "data/resamples/diffexpr"
trace_dir <- "data/resamples/trace_diffexpr"

resample_files <- dir(resample_dir)
trace_files <- dir(trace_dir)

resample_qvals <- sapply(resample_files, function(f) {
  csv <- read_csv(file.path(resample_dir, f))
  csv$q_val
})

trace_qvals <- sapply(trace_files, function(f) {
  csv <- read_csv(file.path(trace_dir, f))
  csv$q_val
})


# Variance analysis -------------------------------------------------------
var_df <- data_frame(DimRed = rowVars(resample_qvals),
                     GPLVM = rowVars(trace_qvals)) %>%
  melt(variable.name = "Variance_source", value.name = "Variance")

ggplot(var_df, aes(x = Variance_source, y = Variance)) + geom_violin()


# Proportion significant analysis -----------------------------------------
prop_df <- data_frame(DimRed = rowMeans(resample_qvals < 0.05),
                      GPLVM = rowMeans(trace_qvals < 0.05))

prop_df %>% 
  melt(variable.name = "Source", value.name = "Proportion_significant") %>%
  ggplot(aes(x = Proportion_significant)) + geom_histogram() +
    facet_wrap(~ Source)

ggplot(prop_df, aes(x = DimRed, y = GPLVM)) + geom_point()


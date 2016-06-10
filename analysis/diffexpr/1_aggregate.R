library(readr)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

study <- args[1]

stopifnot(study %in% c("trapnell", "burns", "shin"))

csv_dir <- file.path("data/diffexpr", study)
csv_files <- dir(csv_dir)

qvals <- lapply(csv_files, function(f) {
  rep <- as.numeric(gsub("de_", "", gsub(".csv", "", f, fixed = TRUE)))
  qv <- read_csv(file.path(csv_dir, f))
  qv <- mutate(qv, rep = rep)
})

qval_df <- do.call(rbind, qvals)

csv_output <- paste0("data/diffexpr/", study, ".csv")
write_csv(qval_df, csv_output)
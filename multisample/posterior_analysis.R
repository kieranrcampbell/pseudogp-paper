
library(rhdf5)
library(reshape2)
library(ggplot2)
library(dplyr)
library(modeest)

get_map <- function(x) mlv(x,method = "HSM")$M

data_dir <- "/net/isi-project/CW010_CAMPBELL_SCNGEN/data/GP/multisample"

h5file <- paste0(data_dir, "/multitrace.h5")

trace <- h5read(h5file, "sample_1/tchain")
trace <- data.frame(trace)
names(trace) <- paste0("t", 1:ncol(trace))
trace$iter <- 1:nrow(trace)

tt <- melt(trace, id.vars = "iter")

ggplot(tt) + geom_line(aes(x = iter, y = value, colour = variable))

burn_period <- as.integer(nrow(trace) / 2)

approx_slopes <- apply(select(trace[burn_period:nrow(trace),], -iter), 2,
                       function(tr) {
                         iter <- (1:length(tr)) / length(tr)
                         fit <- lm(tr ~ iter)
                         coef(fit)[2]
                       })

bad <- which(approx_slopes < -0.1)
qplot(burn_period:nrow(trace), trace[burn_period:nrow(trace), bad], geom = 'line')

h5data <- "/net/isi-scratch/kieran/GP/pseudogp2/data/embeddings.h5"
pst <- h5read(h5data, "monocle/pseudotime")

nsamples <- 50
cors <- sapply(1:nsamples, function(i) {
  sname <- paste0("sample_", i)
  trace <- h5read(h5file, paste0(sname, "/tchain"))
  to_sample <- h5read("/net/isi-scratch/kieran/GP/pseudogp2/data/multisample.h5", 
                      paste0(sname, "/to_sample"))
  map_pseudotime <- apply(trace[burn_period:nrow(trace),], 2, get_map)
  cor(pst[to_sample], map_pseudotime)
})

ggplot(data.frame(x = 1:length(cors), y = cors)) + geom_point(aes(x = x, y = y), size = 3) +
  theme_bw() + scale_y_continuous(breaks = seq(-.9, .9, length.out = 20)) +
  geom_hline(yintercept = c(-.6, .6), color = 'red')



plot_against_pc <- function(i) {
  sname <- paste0("sample_", i)
  trace <- h5read(h5file, paste0(sname, "/tchain"))
  to_sample <- h5read("/net/isi-scratch/kieran/GP/pseudogp2/data/multisample.h5", 
                      paste0(sname, "/to_sample"))
  map_pseudotime <- apply(trace[burn_period:nrow(trace),], 2, get_map)
  qplot(pst[to_sample], map_pseudotime)
}

trouble_points <- abs(cors) < 0.6
trouble_plots <- lapply(which(trouble_points), plot_against_pc)
cowplot::plot_grid(plotlist = plot_list)

boundary_points <- abs(cors) > .6 & abs(cors) < .7
boundary_plots <- lapply(which(boundary_points), plot_against_pc)
cowplot::plot_grid(plotlist = boundary_plots)

ok_points <- abs(cors) > .7 & abs(cors) < .8
ok_plots <- lapply(which(ok_points), plot_against_pc)
cowplot::plot_grid(plotlist = ok_plots)

good_points <- abs(cors) > .8
good_plots <- lapply(which(good_points), plot_against_pc)
cowplot::plot_grid(plotlist = good_plots)




# Testing utility of zero-inflation ---------------------------------------
library(embeddr)


base_dir <- ""
system <- devtools::session_info()$platform$system
if(length(grep("darwin", system)) > 0) {
  # we're on the mac
  base_dir <- "~/mount/"
} else {
  # we're on linux
  base_dir <- "/net/isi-scratch/kieran/"
}


devtools::load_all(paste0(base_dir, "switch/sctools"))
source(paste0(base_dir, "GP/pseudogp2/diffexpr/prep_data.R"))

get_map <- function(x) mlv(x,method = "HSM")$M

data <- load_data(base_dir)
sce <- data$sce ; pst <- data$pst

pst_map <- apply(pst, 2, get_map)
sce$pseudotime <- pst_map

gene <- match("TUBBP1", fData(sce)$gene_short_name)

plot_in_pseudotime(sce[gene,])
x <- exprs(sce)[gene,]

## non-zero inflated model
model <- norm_fit_alt_model(x, pseudotime(sce))

zimodel <- alt_EM_ctrl(x, pseudotime(sce), loglik_tol = 1e-8)

plt1 <- norm_plot_model(model, x, pseudotime(sce))
plt2 <- norm_zi_plot_model(zimodel$par, x, pseudotime(sce))

cowplot::plot_grid(plt1, plt2, ncol = 2)

# Plot the differences ----------------------------------------------------

imputed <- zimodel$x
df <- data.frame(Measured = x, Pseudotime = pseudotime(sce), Imputed = imputed)
dfm <- reshape2::melt(df, id.vars = "Pseudotime", 
                      variable.name = "Type", value.name = "Expression")

dff <- dfm %>% group_by(Pseudotime, Expression) %>% filter(row_number() == 1)

mu_func <- function(t, params) {
  L <- params[1] ; k <- params[2] ; t_0 <- params[3] ; r <- params[4]
  L / (1 + exp(-k*(t - t_0)))
}

blue <- ggthemes::ggthemes_data$fivethirtyeight['blue']
red <- ggthemes::ggthemes_data$fivethirtyeight['red']

ggplot(dff) + geom_point(aes(x = Pseudotime, y = Expression, colour = Type), size =2.5) +
  ggthemes::scale_colour_fivethirtyeight() +
  stat_function(fun = mu_func, args = list(model$par), color=blue) +
  stat_function(fun = mu_func, args = list(zimodel$par), color = red) 
  
ggsave(filename = paste0(base_dir, "GP/pseudogp2/diffexpr/zi.png"),
       width = 3.2, height = 1.8, scale = 3)


# Quick play with flexmix -------------------------------------------------
library(flexmix)
pst <- pseudotime(sce)
plot(pst, x)
df <- data.frame(t = pst, x = x)
df$xp <- df$x + replicate(length(df$x), jitter(0))

cc <- FLXPmultinom(~ t)
m <- flexmix(xp ~ I(1 / (1 + exp(-t))),  data = df, k = 2, concomitant = cc)

df$cluster <- clusters(m)
ggplot(df, aes(x = t, y = x, colour = cluster)) + geom_point()






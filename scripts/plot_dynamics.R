# Figure 1: short trajectories of the Zig-Zag, Bouncy Particle, and Boomerang
# samplers on a two-dimensional banana-shaped target.
#
# Run from the project root:
#     Rscript scripts/plot_dynamics.R

library(ggplot2)
library(patchwork)
library(PDMPSamplersR)
source("scripts/figure_theme.R")   # theme_paper(), save_figure()

# ---- Target -------------------------------------------------------------------

# Banana density: x ~ Normal(mu, 1/(2a)), y | x ~ Normal(x^2, 1/(2b)).
mu <- 0
a <- 1 / 2
b <- 2
dbanana <- function(x, y) dnorm(x, mu, 1 / sqrt(2 * a)) * dnorm(y, x^2, 1 / sqrt(2 * b))
# Gradient of the negative log density, as the samplers expect.
grad_neg_log_banana <- function(theta) {
  r <- theta[2] - theta[1]^2
  -c(-2 * a * (theta[1] - mu) + 4 * b * theta[1] * r, -2 * b * r)
}

grid <- expand.grid(x = seq(-3, 3, length.out = 256), y = seq(-1.5, 4.5, length.out = 256))
grid$z <- dbanana(grid$x, grid$y)

# ---- Trajectories ---------------------------------------------------------------

x0 <- c(-1, 1)
flows <- c("Zig-Zag" = "ZigZag", "Bouncy Particle" = "BouncyParticle",
           "Boomerang" = "AdaptiveBoomerang")
dt <- 0.01
n_shown <- 5000                       # the first 50 units of PDMP time
paths <- do.call(rbind, lapply(names(flows), function(name) {
  trace <- pdmp_sample(grad_neg_log_banana, 2, flow = flows[[name]],
                       algorithm = "GridThinningStrategy", x0 = x0, T = 5000,
                       seed = 42, show_progress = FALSE)
  s <- discretize(trace, dt)[seq_len(n_shown), ]
  data.frame(sampler = name, x = s[, 1], y = s[, 2])
}))
paths$sampler <- factor(paths$sampler, levels = names(flows))

# ---- Figure ---------------------------------------------------------------------

p <- ggplot(grid, aes(x, y)) +
  geom_contour(aes(z = z), linewidth = 0.25, bins = 15, colour = "grey55") +
  geom_path(data = paths, colour = "#08519c", linewidth = 0.35) +
  geom_point(data = paths[!duplicated(paths$sampler), ], shape = 21, size = 2,
             fill = "grey80", colour = "black") +
  facet_wrap(~sampler, nrow = 1) +
  scale_x_continuous(breaks = seq(-2, 2, 2), labels = NULL) +
  scale_y_continuous(breaks = seq(-2, 4, 2), labels = NULL) +
  coord_cartesian(xlim = c(-3, 3), ylim = c(-1.5, 4.5), expand = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_paper() +
  theme(aspect.ratio = 1, axis.ticks.length = unit(1, "mm"))

save_figure(p, "dynamics", height_cm = 6)

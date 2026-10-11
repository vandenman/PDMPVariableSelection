# Figure of the empirical RI-CLPM example (scripts/empirical_clpm.R).
#
# Run from the project root after the analysis:
#     Rscript scripts/plot_empirical_clpm.R
# Panel A: the cross-lagged network of the Boomerang fits under the dependent
# prior (both chains): every cross-lag with posterior inclusion probability
# above 0.5, drawn with width proportional to its marginal posterior mean.
# Panels B and C: the trajectories of two cross-lags with a common outcome and
# strongly correlated predictors, from the shorter run that keeps its trajectory.

library(qgraph)
library(grid)
source("scripts/figure_theme.R")   # theme_paper(), save_figure()

out_dir <- "results/empirical_liss_health"
labels <- c(anxious = "anxious", down = "down", calm_r = "not\ncalm",
            depressed = "depressed", happy_r = "not\nhappy", poor_health = "poor\nhealth",
            hinder_daily = "hinders\ndaily", hinder_social = "hinders\nsocial",
            pain = "pain", fatigue = "fatigue", sleep = "sleep\nproblems")
mood <- c("anxious", "down", "calm_r", "depressed", "happy_r")
pair <- c("depressed", "down")      # predictors of the shown pair
pair_outcome <- "fatigue"

# ---- Network ------------------------------------------------------------------

source("scripts/simulation_common.R")   # fit_parameters()
variables <- names(labels)
cross <- do.call(rbind, lapply(1:2, function(chain) {
  f <- readRDS(file.path(out_dir, "fits", sprintf("dependent_AdaptiveBoomerang_chain%d.rds", chain)))
  p <- fit_parameters(f)
  ij <- regmatches(p$parameter, regexec("^A\\.([0-9]+)\\.([0-9]+)$", p$parameter))
  keep <- lengths(ij) == 3
  data.frame(outcome = variables[as.integer(sapply(ij[keep], `[`, 2))],
             predictor = variables[as.integer(sapply(ij[keep], `[`, 3))],
             pip = p$pip[keep], mean = p$mean[keep])
}))
cross <- cross[cross$outcome != cross$predictor, ]
edges <- aggregate(cbind(pip, mean) ~ predictor + outcome, cross, mean)
W <- matrix(0, length(variables), length(variables), dimnames = list(variables, variables))
shown <- edges[edges$pip > 0.5, ]
W[cbind(shown$predictor, shown$outcome)] <- shown$mean   # row = predictor, column = outcome

# ---- Trajectories ---------------------------------------------------------------

trajectory <- readRDS(file.path(out_dir, "trajectory.rds"))
column <- function(predictor, outcome)
  sprintf("A.%d.%d", match(outcome, trajectory$variables), match(predictor, trajectory$variables))
time <- trajectory$trajectory[, "time"] - trajectory$trajectory[1, "time"]
paths <- sapply(pair, function(p) trajectory$trajectory[, column(p, pair_outcome)])
# Show the 150 time units with the most swaps between the two cross-lags.
swaps <- cumsum(c(0, abs(diff(paths[, 1] != 0))))
width <- 150
starts <- seq(0, max(time) - width, by = 5)
best <- starts[which.max(sapply(starts, function(s)
  diff(range(swaps[time >= s & time <= s + width]))))]
window <- time >= best & time <= best + width

# ---- Figure -----------------------------------------------------------------------

# Panels B and C (ggplot2): the trajectories of the two cross-lags in the window.
trace_panel <- function(k, tag) {
  d <- data.frame(time = time[window] - best, value = paths[window, k])
  ggplot(d, aes(time, value)) +
    geom_hline(yintercept = 0, colour = "grey60", linetype = 2) +
    geom_line(colour = "#08519c", linewidth = 0.35) +
    scale_y_continuous(limits = range(c(0, paths[window, ]))) +
    labs(x = if (k == length(pair)) "PDMP time" else NULL,
         y = bquote(.(pair[k]) %->% .(pair_outcome)), tag = tag) +
    theme_paper()
}

# Panel A (qgraph) is drawn with base graphics in the left part of the page, and
# the ggplot2 panels are printed into grid viewports on the right.
draw <- function() {
  grid.newpage()
  suppressWarnings(par(fig = c(0, 0.5, 0, 1), mar = c(0, 0, 0, 0), new = TRUE))
  qgraph(W, layout = "circle", directed = TRUE, labels = labels[variables],
         groups = list(Mood = match(mood, variables),
                       `Physical health` = match(setdiff(variables, mood), variables)),
         color = c("#9ecae1", "#fdd0a2"), legend = FALSE, label.scale = FALSE,
         label.cex = 0.6, vsize = 14, maximum = max(abs(W)), fade = FALSE,
         posCol = "#08519c", negCol = "#a50f15", asize = 4, edge.width = 1.5,
         mar = c(3, 3, 3, 3))
  grid.text("A", x = unit(2, "mm"), y = unit(1, "npc") - unit(3, "mm"),
            just = c("left", "top"), gp = gpar(fontsize = 12))
  print(trace_panel(1, "B"), vp = viewport(x = 0.75, y = 0.75, width = 0.5, height = 0.5))
  print(trace_panel(2, "C"), vp = viewport(x = 0.75, y = 0.25, width = 0.5, height = 0.5))
}
size <- c(width = figure_width_cm, height = 8.5) / 2.54      # inches
cairo_pdf(file.path(figure_dir, "empirical_liss_health.pdf"), width = size[[1]], height = size[[2]])
draw(); invisible(dev.off())
png(file.path(figure_dir, "empirical_liss_health.png"), width = size[[1]], height = size[[2]],
    units = "in", res = 300)
draw(); invisible(dev.off())
cat(sprintf("%d cross-lags shown (PIP > 0.5); trajectory window %.0f-%.0f with %d swaps\n",
            nrow(shown), best, best + width, diff(range(swaps[window]))))


# ---- Bivariate trajectory ---------------------------------------------------------

# The two cross-lags against each other over part of the window (time 25 to 75
# of panels B and C), to show how the process moves between the two explanations.
part <- window & time - best >= 25 & time - best <= 75
d <- data.frame(time = time[part] - best, x = paths[part, 1], y = paths[part, 2])
p <- ggplot(d, aes(x, y, colour = time)) +
  geom_hline(yintercept = 0, colour = "grey60", linetype = 2) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
  geom_path(linewidth = 0.4) +
  geom_point(data = d[c(1, nrow(d)), ], shape = c(1, 4), size = 2.5, colour = "black") +
  scale_colour_viridis_c(option = "D", end = 0.9) +
  coord_equal() +
  labs(x = bquote(.(pair[1]) %->% .(pair_outcome)), y = bquote(.(pair[2]) %->% .(pair_outcome)),
       colour = "PDMP time") +
  theme_paper() + theme(legend.position = "right")
save_figure(p, "empirical_liss_health_bivariate", height_cm = 9, width_cm = 12)

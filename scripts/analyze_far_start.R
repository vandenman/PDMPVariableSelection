# Time to the typical set from a start far from the posterior mode.
#
# Run from the project root after scripts/simulation_far_start.R:
#     Rscript scripts/analyze_far_start.R
#
# For every batch of the PDMP warmup or the NIMBLE burn-in, the batch mean of
# each edge is compared with the reference posterior of its data set, the
# pooled estimates of the mode-start fits of the OMRF study (all methods):
#   distance = mean over the edges of |batch mean - posterior mean| / posterior sd.
# The time to the typical set is the first time, in seconds of computation, at
# which the distance falls below 1 (time interpolated linearly within the
# warmup or burn-in). Writes
#   simulation_results/far_start/far_start_typical_set.csv   per fit
#   ../proposal/figures/far_start_typical_set.pdf

source("scripts/simulation_common.R")
source("scripts/figure_theme.R")   # theme_paper(), save_figure(), sampler_labels

dataset <- function(task) sprintf("N%d_P%d_%s_rep%02d__%s", task$N, task$P,
                                  task$structure, task$rep, task$prior)

# Reference posterior per data set from the mode-start study fits.
reference <- list()
for (f in list.files("simulation_results/omrf_study/fits", full.names = TRUE,
                     pattern = "^N200_P20_random_rep[0-9]+__independent_")) {
  r <- readRDS(f)
  if (!is.null(r$error)) next
  p <- r$parameters[grepl("^interactions_0", r$parameters$parameter), ]
  reference[[dataset(r$task)]] <- rbind(reference[[dataset(r$task)]],
    data.frame(key = gsub("\\W", "", p$parameter), mean = p$mean, sd = p$sd))
}
reference <- lapply(reference, function(r) data.frame(
  mean = tapply(r$mean, r$key, mean), sd = tapply(r$sd, r$key, median, na.rm = TRUE)))

per_fit <- do.call(rbind, lapply(list.files("simulation_results/far_start/fits",
                                            full.names = TRUE), function(f) {
  r <- readRDS(f)
  row <- data.frame(distance = r$task$start_distance, rep = r$task$rep,
                    sampler = r$task$sampler, started = is.null(r$error),
                    time_typical = NA_real_, warmup_seconds = NA_real_)
  if (is.null(r$error)) {
    w <- r$convergence$warmup
    ref <- reference[[dataset(r$task)]]
    k <- match(gsub("\\W", "", rownames(w$x)), rownames(ref))
    distance <- colMeans(abs(sweep(w$x, 2L, w$weight, "/") - ref$mean[k]) / ref$sd[k])
    time <- w$seconds * cumsum(w$weight) / sum(w$weight)
    row$time_typical <- time[which(distance < 1)[1]]
    row$warmup_seconds <- w$seconds
  }
  row
}))
write.csv(per_fit, "simulation_results/far_start/far_start_typical_set.csv", row.names = FALSE)

per_fit$method <- factor(sampler_labels[per_fit$sampler], sampler_labels)
summary <- do.call(rbind, lapply(split(per_fit, list(per_fit$distance, per_fit$method),
                                       drop = TRUE), function(g) data.frame(
  distance = g$distance[1], method = g$method[1], fits = nrow(g),
  started = sum(g$started), reached = sum(!is.na(g$time_typical)),
  median = median(g$time_typical, na.rm = TRUE),
  min = suppressWarnings(min(g$time_typical, na.rm = TRUE)),
  max = suppressWarnings(max(g$time_typical, na.rm = TRUE)))))
print(summary, digits = 3, row.names = FALSE)

# Figure: time to the typical set per fit, with the median; fits that did not
# start (NIMBLE) or did not reach the typical set within the warmup are counted
# above the points.
d <- per_fit[!is.na(per_fit$time_typical), ]
counts <- aggregate(cbind(not_started = !started, not_reached = started & is.na(time_typical)) ~
                      distance + method, per_fit, sum)
counts$label <- ifelse(counts$not_started > 0,
                       paste(counts$not_started, "of 10\ndid not start"),
                       paste(counts$not_reached, "of 10\nnot reached"))
counts <- counts[counts$not_started + counts$not_reached > 0, ]
panels <- paste("Distance", sort(unique(per_fit$distance)))
d$panel <- factor(paste("Distance", d$distance), panels)
counts$panel <- factor(paste("Distance", counts$distance), panels)
p <- ggplot(d, aes(method, time_typical)) +
  geom_point(position = position_jitter(width = 0.12, height = 0, seed = 1),
             size = 1.2, alpha = 0.6) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, linewidth = 0.3) +
  geom_text(data = counts, aes(method, Inf, label = label),
            vjust = 1.3, size = 2.6, lineheight = 0.9) +
  facet_wrap(~ panel) +
  scale_x_discrete(drop = FALSE, labels = function(l) sub(" ", "\n", l)) +
  blank_125(d, "time_typical", "panel", iqr = FALSE) +
  scale_y_log10(breaks = breaks_125, labels = labels_125, expand = expansion(mult = c(0.05, 0.35))) +
  labs(x = NULL, y = "Time to the typical set (seconds)") +
  theme_paper()
save_figure(p, "far_start_typical_set", height_cm = 7)

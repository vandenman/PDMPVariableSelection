# Analysis of the OMRF simulation study (scripts/simulation_omrf.R).
#
# Run from the project root after the study has finished:
#     Rscript scripts/analyze_omrf_study.R
# It reads simulation_results/omrf_study/omrf_{fits,parameters}.csv and writes
#   summary tables  simulation_results/omrf_study/summaries/
#   figures         simulation_results/omrf/final/plots/ (the manuscript path)
# The sections are: 1. per-fit measures  2. summary tables  3. figures.

library(ggplot2)

study_dir <- "simulation_results/omrf_study"
summary_dir <- file.path(study_dir, "summaries")
plot_dir <- "simulation_results/omrf/final/plots"
dir.create(summary_dir, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

fits <- read.csv(file.path(study_dir, "omrf_fits.csv"))
parameters <- read.csv(file.path(study_dir, "omrf_parameters.csv"))
stopifnot(all(is.na(fits$error)))

sampler_labels <- c(AdaptiveBoomerang = "Boomerang", PreconditionedBPS = "BPS",
                    PreconditionedZigZag = "Zig-Zag", nimble = "NIMBLE RJ")
prior_labels <- c(independent = "Independent prior", node = "Node prior")
# One method per sampler, gradient and prior.
fits$method <- with(fits, ifelse(sampler == "nimble",
  paste0("NIMBLE RJ, ", prior),
  paste0(sampler_labels[sampler], ", ", gradient,
         ifelse(prior == "node", ", node prior", ""))))

# ---- 1. Per-fit measures ------------------------------------------------------

# Edge interactions only (NIMBLE also reports its node scales).
edges <- parameters[grepl("^interactions_", parameters$parameter), ]
edges <- merge(edges[, c("id", "parameter", "mean", "pip", "ess_mean", "ess_pip",
                         "true_value")],
               fits[, c("id", "P", "N", "structure", "rep", "prior", "sampler",
                        "gradient", "method")],
               by = "id")
edges$active <- edges$true_value != 0

# Area under the ROC curve of the PIPs (Mann-Whitney form, ties count 1/2).
auc <- function(pip, active) {
  if (all(active) || !any(active)) return(NA_real_)
  mean(outer(pip[active], pip[!active], ">") +
         0.5 * outer(pip[active], pip[!active], "=="))
}

per_fit <- do.call(rbind, lapply(split(edges, edges$id), function(d) {
  selected <- d$pip >= 0.5
  data.frame(
    id = d$id[1],
    median_ess_mean = median(d$ess_mean, na.rm = TRUE),
    median_ess_pip = median(d$ess_pip, na.rm = TRUE),   # undefined when PIP = 1
    rmse = sqrt(mean((d$mean - d$true_value)^2)),
    auc = auc(d$pip, d$active),
    sensitivity = mean(selected[d$active]),
    specificity = mean(!selected[!d$active]),
    model_size = sum(d$pip),           # posterior expected number of edges
    true_size = sum(d$active))
}))
per_fit <- merge(fits, per_fit, by = "id")
# Efficiency: effective sample size of the whole retained trace (PDMP: retained
# trajectory after warmup; NIMBLE: retained draws) per second of wall time after
# compilation. PDMP: `warmup` + `sampling` (the sampler's phase timers; with
# subsampling also the full-gradient run that supplies the anchors). NIMBLE:
# `sampling` is runMCMC, burn-in plus retained iterations. Every fit ran on one
# core.
nimble <- per_fit$sampler == "nimble"
per_fit$run_time <- ifelse(nimble, per_fit$sampling, per_fit$warmup + per_fit$sampling)
per_fit$ess_per_second <- per_fit$median_ess_mean / per_fit$run_time
per_fit$pip_ess_per_second <- per_fit$median_ess_pip / per_fit$run_time

# Agreement with the NIMBLE reference fit of the same data set and prior (P = 10).
reference <- edges[edges$sampler == "nimble",
                   c("P", "N", "structure", "rep", "prior", "parameter", "mean", "pip")]
names(reference)[7:8] <- c("ref_mean", "ref_pip")
agreement <- merge(edges[edges$sampler != "nimble", ], reference,
                   by = c("P", "N", "structure", "rep", "prior", "parameter"))
agreement_fit <- aggregate(
  cbind(pip_diff = abs(pip - ref_pip), mean_diff = abs(mean - ref_mean)) ~
    id + P + N + structure + prior + method, agreement, mean)

# ---- 2. Summary tables --------------------------------------------------------

med <- function(formula, data) aggregate(formula, data, median, na.rm = TRUE,
                                          na.action = na.pass)

efficiency <- med(cbind(ess_per_second, pip_ess_per_second, run_time, events,
                        boundary_recoveries) ~ method + P + N, per_fit)
selection <- med(cbind(auc, sensitivity, specificity, rmse, model_size, true_size) ~
                   method + P + N, per_fit)
agree <- med(cbind(pip_diff, mean_diff) ~ method + N, agreement_fit)

write.csv(per_fit, file.path(summary_dir, "omrf_per_fit.csv"), row.names = FALSE)
write.csv(efficiency, file.path(summary_dir, "omrf_efficiency.csv"), row.names = FALSE)
write.csv(selection, file.path(summary_dir, "omrf_selection.csv"), row.names = FALSE)
write.csv(agree, file.path(summary_dir, "omrf_agreement_with_nimble.csv"),
          row.names = FALSE)

print(efficiency, digits = 3)
print(selection, digits = 3)
print(agree, digits = 3)

# ---- 3. Figures -----------------------------------------------------------------

label_fits <- function(d) {
  d$sampler_label <- factor(sampler_labels[d$sampler], levels = sampler_labels)
  d$prior_label <- factor(prior_labels[d$prior], levels = prior_labels)
  d$gradient_label <- factor(ifelse(is.na(d$gradient), "MCMC", d$gradient),
                             levels = c("full", "subsampled", "MCMC"))
  d$N_label <- factor(paste0("N = ", d$N))
  d$P_label <- factor(paste0("p = ", d$P), levels = paste0("p = ", sort(unique(d$P))))
  d
}
gradient_colours <- c(full = "#0072B2", subsampled = "#E69F00", MCMC = "#009E73")
theme_set(theme_bw(base_size = 11) +
            theme(panel.grid.minor = element_blank(),
                  axis.text.x = element_text(angle = 25, hjust = 1)))
save_plot <- function(plot, name, width = 9, height = 6) {
  ggsave(file.path(plot_dir, paste0(name, ".pdf")), plot, width = width,
         height = height)
  ggsave(file.path(plot_dir, paste0(name, ".png")), plot, width = width,
         height = height, dpi = 300)
}

# Median over replications and structures, one point per sampler and gradient.
condition_plot <- function(y, ylab, log = FALSE) {
  d <- label_fits(per_fit)
  p <- ggplot(d, aes(sampler_label, .data[[y]], colour = gradient_label)) +
    stat_summary(fun = median, geom = "point", size = 2.4,
                 position = position_dodge(width = 0.5)) +
    facet_grid(P_label ~ prior_label + N_label, scales = "free_y") +
    scale_colour_manual(values = gradient_colours) +
    labs(x = NULL, y = ylab, colour = "Gradient")
  if (log) p + scale_y_log10() else p
}

save_plot(condition_plot("ess_per_second", "Median edge ESS per second", log = TRUE),
          "omrf_ess_per_second_by_sampler")
save_plot(condition_plot("run_time", "Warmup + sampling time (seconds)", log = TRUE),
          "omrf_runtime_by_sampler")
save_plot(condition_plot("rmse", "RMSE over edge interactions"), "omrf_accuracy_rmse")
save_plot(condition_plot("auc", "AUC of the inclusion probabilities"), "omrf_auc")

# PIPs of each PDMP method against the NIMBLE reference (P = 10).
d <- label_fits(agreement)
p <- ggplot(d, aes(ref_pip, pip)) +
  geom_abline(colour = "grey60") +
  geom_point(alpha = 0.3, size = 0.8, colour = "#0072B2") +
  facet_wrap(~ method, nrow = 2) +
  coord_equal() +
  labs(x = "PIP, NIMBLE reference", y = "PIP, PDMP sampler")
save_plot(p, "omrf_pip_agreement_with_nimble", width = 10, height = 5.5)

message("Summaries: ", summary_dir, "\nFigures: ", plot_dir)

# The manuscript figure, at print size: one point per sampler and gradient,
# the median over replications and graph structures.
source("scripts/figure_theme.R")   # theme_paper(), save_figure(), sampler_labels
slab_labels <- c(independent = "Independent slab", node = "Node-scale slab")
method_labels <- c(full = "Full gradient", subsampled = "Subsampled gradient",
                   MCMC = "NIMBLE")
d <- per_fit
d$sampler_label <- factor(sampler_labels[d$sampler], levels = sampler_labels)
d$gradient_label <- factor(method_labels[ifelse(d$sampler == "nimble", "MCMC", d$gradient)],
                           levels = method_labels)
d$column_label <- factor(paste0(slab_labels[d$prior], "\nN = ", d$N),
                         levels = c(outer(slab_labels, paste0("\nN = ", sort(unique(d$N))), paste0)))
d$P_label <- factor(paste0("p = ", d$P), levels = paste0("p = ", sort(unique(d$P))))
p <- ggplot(d, aes(sampler_label, ess_per_second, colour = gradient_label,
                   shape = gradient_label)) +
  stat_summary(fun.data = median_iqr, geom = "linerange", linewidth = 0.4,
               position = position_dodge(width = 0.4)) +
  stat_summary(fun = median, geom = "point", size = 1.8,
               position = position_dodge(width = 0.4)) +
  facet_grid(P_label ~ column_label, scales = "free", space = "free_x") +
  blank_125(d, "ess_per_second", c("P_label", "column_label"), c("sampler_label", "gradient_label")) +
  scale_y_log10(breaks = breaks_125, labels = labels_125) +
  scale_colour_manual(values = c("#0072B2", "#E69F00", "#009E73")) +
  scale_shape_manual(values = c(16, 17, 15)) +
  labs(x = NULL, y = "Median ESS of the edges per second", colour = NULL, shape = NULL) +
  theme_paper() + theme(axis.text.x = element_text(angle = 40, hjust = 1))
save_figure(p, "omrf_ess_per_second_by_sampler", height_cm = 15)

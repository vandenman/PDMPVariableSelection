# Time to convergence in the simulation studies (Vogels et al., 2024).
#
# Run from the project root after a study has finished:
#     Rscript scripts/analyze_convergence.R clpm   (or omrf)
# For every fit, the running inclusion probabilities and model-averaged
# posterior means over computation time come from the stored batch integrals
# (PDMP) or batch sums of the draws (NIMBLE), from the start of the run
# (convergence_trace()). The time to convergence is the first time at which
#   the AUC of the inclusion probabilities is within 0.01 of its final value,
#   the RMSE of the posterior means is within 5% of its final value
# (convergence_times()). Writes
#   simulation_results/<study>_study/summaries/<study>_convergence.csv  per fit
#   simulation_results/<study>_study/summaries/<study>_convergence_summary.csv
#   ../proposal/figures/<study>_time_to_convergence.pdf

study <- commandArgs(TRUE)[1]
source("scripts/simulation_common.R")   # convergence_times(), fit_parameters()
source("scripts/figure_theme.R")        # theme_paper(), save_figure(), sampler_labels

study_dir <- sprintf("simulation_results/%s_study", study)
summary_dir <- file.path(study_dir, "summaries")
dir.create(summary_dir, showWarnings = FALSE)

per_fit <- do.call(rbind, lapply(list.files(file.path(study_dir, "fits"),
                                            full.names = TRUE), function(file) {
  f <- readRDS(file)
  if (!is.null(f$error) || is.null(f$convergence)) return(NULL)
  p <- f$parameters
  p <- p[if (study == "clpm") p$selectable & !is.na(p$true_value) else
           grepl("^interactions_0", p$parameter), ]
  times <- convergence_times(f$convergence, setNames(p$true_value, p$parameter))
  cbind(as.data.frame(f$task[intersect(names(f$task),
          c("id", "N", "P", "pattern", "structure", "rep", "prior", "sampler", "gradient"))]),
        as.list(times))
}))
write.csv(per_fit, file.path(summary_dir, sprintf("%s_convergence.csv", study)),
          row.names = FALSE)

per_fit$gradient <- per_fit$gradient %||% NA
per_fit$method <- with(per_fit, ifelse(is.na(gradient), sampler_labels[sampler],
  paste0(sampler_labels[sampler], ifelse(gradient == "subsampled", " (subsampled)", ""))))
by <- intersect(c("method", "prior", "P", "N"), names(per_fit))
summary <- aggregate(cbind(time_auc, time_rmse, total_time) ~ .,
                     per_fit[, c(by, "time_auc", "time_rmse", "total_time")],
                     median, na.rm = TRUE, na.action = na.pass)  # no AUC in null conditions
write.csv(summary, file.path(summary_dir, sprintf("%s_convergence_summary.csv", study)),
          row.names = FALSE)
print(summary, digits = 3, row.names = FALSE)

# Figure: time to convergence of the AUC per sampler and condition (median and
# interquartile range over the fits), laid out as the efficiency figures.
prior_labels <- if (study == "omrf")
  c(independent = "Independent slab", node = "Node-scale slab") else
  c(independent = "Independent", dependent = "Dependent")
method_labels <- c(full = "Full gradient", subsampled = "Subsampled gradient",
                   MCMC = "NIMBLE")
d <- per_fit[!is.na(per_fit$time_auc), ]
d$sampler_label <- factor(sampler_labels[d$sampler], levels = sampler_labels)
d$gradient_label <- factor(method_labels[ifelse(d$sampler == "nimble", "MCMC",
                                                ifelse(is.na(d$gradient), "full", d$gradient))],
                           levels = method_labels)
d$column_label <- factor(paste0(prior_labels[d$prior], "\nN = ", d$N),
  levels = c(outer(prior_labels, paste0("\nN = ", sort(unique(d$N))), paste0)))
d$P_label <- factor(paste0("p = ", d$P), levels = paste0("p = ", sort(unique(d$P))))
p <- ggplot(d, aes(sampler_label, time_auc, colour = gradient_label, shape = gradient_label)) +
  stat_summary(fun.data = median_iqr, geom = "linerange", linewidth = 0.4,
               position = position_dodge(width = 0.4)) +
  stat_summary(fun = median, geom = "point", size = 1.8,
               position = position_dodge(width = 0.4)) +
  facet_grid(P_label ~ column_label, scales = "free", space = "free_x") +
  blank_125(d, "time_auc", c("P_label", "column_label"), c("sampler_label", "gradient_label")) +
  scale_y_log10(breaks = breaks_125, labels = labels_125) +
  scale_colour_manual(values = c("#0072B2", "#E69F00", "#009E73"), drop = TRUE) +
  scale_shape_manual(values = c(16, 17, 15), drop = TRUE) +
  labs(x = NULL, y = "Time to AUC convergence (seconds)", colour = NULL, shape = NULL) +
  theme_paper() + theme(axis.text.x = element_text(angle = 40, hjust = 1))
save_figure(p, sprintf("%s_time_to_convergence", study), height_cm = 15)


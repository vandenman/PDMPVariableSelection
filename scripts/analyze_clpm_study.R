# Analysis of the RI-CLPM simulation study (scripts/simulation_clpm.R).
#
# Run from the project root after the study has finished:
#     Rscript scripts/analyze_clpm_study.R
# It reads simulation_results/clpm_study/clpm_{fits,parameters}.csv and writes
#   summary tables  simulation_results/clpm_study/summaries/
#   figures         simulation_results/clpm/final/plots/ (the manuscript path)
# The sections are: 1. per-fit measures  2. summary tables  3. figures.

library(ggplot2)

study_dir <- "simulation_results/clpm_study"
summary_dir <- file.path(study_dir, "summaries")
plot_dir <- "simulation_results/clpm/final/plots"
dir.create(summary_dir, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

fits <- read.csv(file.path(study_dir, "clpm_fits.csv"))
parameters <- read.csv(file.path(study_dir, "clpm_parameters.csv"))
stopifnot(all(is.na(fits$error)))

source("scripts/figure_theme.R")   # theme_paper(), sampler_labels, save_figure()
pattern_labels <- c(null = "Null", sender = "Sender", mixed = "Mixed")
prior_labels <- c(independent = "Independent prior", dependent = "Dependent prior")

# ---- 1. Per-fit measures ------------------------------------------------------

# Cross-lagged effects only, keyed by (outcome, predictor). The PDMP output
# names them A.i.j, NIMBLE A[i,j].
cross <- parameters[parameters$selectable, ]
cross$key <- gsub("^A[.\\[]([0-9]+)[.,]([0-9]+)\\]?$", "\\1,\\2", cross$parameter)
cross <- merge(cross[, c("id", "key", "mean", "pip", "ess_mean", "ess_pip",
                         "true_value")],
               fits[, c("id", "P", "N", "pattern", "rep", "prior", "sampler")],
               by = "id")
cross$active <- cross$true_value != 0

# Area under the ROC curve of the PIPs (Mann-Whitney form, ties count 1/2).
auc <- function(pip, active) {
  if (all(active) || !any(active)) return(NA_real_)
  mean(outer(pip[active], pip[!active], ">") +
         0.5 * outer(pip[active], pip[!active], "=="))
}

per_fit <- do.call(rbind, lapply(split(cross, cross$id), function(d) {
  selected <- d$pip >= 0.5
  data.frame(
    id = d$id[1],
    # ESS is undefined for an effect that never left zero in the retained run.
    median_ess_mean = median(d$ess_mean, na.rm = TRUE),
    never_free = sum(is.na(d$ess_mean)),
    median_ess_pip = median(d$ess_pip, na.rm = TRUE),   # undefined when PIP = 1
    rmse = sqrt(mean((d$mean - d$true_value)^2)),
    auc = auc(d$pip, d$active),
    sensitivity = if (any(d$active)) mean(selected[d$active]) else NA,
    specificity = mean(!selected[!d$active]),
    false_positives = sum(selected & !d$active),
    null_decisions = sum(!d$active),
    model_size = sum(d$pip))           # posterior expected number of cross-lags
}))
per_fit <- merge(fits, per_fit, by = "id")
# Efficiency: effective sample size of the whole retained trace (PDMP: retained
# trajectory after warmup; NIMBLE: retained draws) per second of wall time after
# compilation. PDMP: `warmup` + `sampling` (the sampler's phase timers). NIMBLE:
# `sampling` is runMCMC, burn-in plus retained iterations. Every fit ran on one
# core.
nimble <- per_fit$sampler == "nimble"
per_fit$run_time <- ifelse(nimble, per_fit$sampling, per_fit$warmup + per_fit$sampling)
per_fit$ess_per_second <- per_fit$median_ess_mean / per_fit$run_time
per_fit$pip_ess_per_second <- per_fit$median_ess_pip / per_fit$run_time

# Agreement with the NIMBLE reference fit of the same data set and prior.
reference <- cross[cross$sampler == "nimble",
                   c("P", "N", "pattern", "rep", "prior", "key", "mean", "pip")]
names(reference)[7:8] <- c("ref_mean", "ref_pip")
agreement <- merge(cross[cross$sampler != "nimble", ], reference,
                   by = c("P", "N", "pattern", "rep", "prior", "key"))
agreement_fit <- aggregate(
  cbind(pip_diff = abs(pip - ref_pip), mean_diff = abs(mean - ref_mean)) ~
    id + P + N + pattern + prior + sampler, agreement, mean)

# ---- 2. Summary tables --------------------------------------------------------

med <- function(formula, data) aggregate(formula, data, median, na.rm = TRUE,
                                          na.action = na.pass)

efficiency <- med(cbind(ess_per_second, pip_ess_per_second, run_time, total,
                        events, gradient_calls) ~ sampler + P + N, per_fit)
selection <- med(cbind(auc, sensitivity, specificity, rmse) ~ sampler + P + N,
                 per_fit[per_fit$pattern != "null", ])
false_pos <- aggregate(cbind(false_positives, null_decisions) ~ sampler + P + N,
                       per_fit, sum)
null_size <- med(model_size ~ prior + P + N + sampler, per_fit[per_fit$pattern == "null", ])
agree <- med(cbind(pip_diff, mean_diff) ~ sampler + prior + P + N, agreement_fit)

write.csv(per_fit, file.path(summary_dir, "clpm_per_fit.csv"), row.names = FALSE)
write.csv(efficiency, file.path(summary_dir, "clpm_efficiency.csv"), row.names = FALSE)
write.csv(selection, file.path(summary_dir, "clpm_selection.csv"), row.names = FALSE)
write.csv(false_pos, file.path(summary_dir, "clpm_false_positives.csv"), row.names = FALSE)
write.csv(null_size, file.path(summary_dir, "clpm_null_model_size.csv"), row.names = FALSE)
write.csv(agree, file.path(summary_dir, "clpm_agreement_with_nimble.csv"),
          row.names = FALSE)

print(efficiency, digits = 3)
print(selection, digits = 3)
print(false_pos)
print(null_size, digits = 3)
print(agree, digits = 3)

# ---- 3. Figures -----------------------------------------------------------------

label_fits <- function(d) {
  d$sampler_label <- factor(sampler_labels[d$sampler], levels = sampler_labels)
  d$pattern_label <- factor(pattern_labels[d$pattern], levels = pattern_labels)
  d$prior_label <- factor(prior_labels[d$prior], levels = prior_labels)
  d$N_label <- factor(paste0("N = ", d$N))
  d$column_label <- factor(paste0(sub(" prior", "", d$prior_label), ", N = ", d$N),
                           levels = c(outer(sub(" prior", "", prior_labels), paste0(", N = ", sort(unique(d$N))), paste0)))
  d$P_label <- factor(paste0("p = ", d$P), levels = paste0("p = ", sort(unique(d$P))))
  d
}
pattern_colours <- c(Null = "#0072B2", Sender = "#E69F00", Mixed = "#009E73")
pattern_shapes <- c(Null = 16, Sender = 17, Mixed = 15)
theme_set(theme_paper() + theme(axis.text.x = element_text(angle = 40, hjust = 1)))
save_plot <- function(plot, name, width = 9, height = 6) {
  ggsave(file.path(plot_dir, paste0(name, ".pdf")), plot, width = width,
         height = height)
  ggsave(file.path(plot_dir, paste0(name, ".png")), plot, width = width,
         height = height, dpi = 300)
}

# Median over replications, one point per sampler, pattern, prior and N.
condition_plot <- function(y, ylab, log = FALSE) {
  d <- label_fits(per_fit)
  p <- ggplot(d, aes(sampler_label, .data[[y]], colour = pattern_label,
                     shape = pattern_label)) +
    stat_summary(fun.data = median_iqr, geom = "linerange", linewidth = 0.4,
                 position = position_dodge(width = 0.4)) +
    stat_summary(fun = median, geom = "point", size = 1.8,
                 position = position_dodge(width = 0.4)) +
    facet_grid(P_label ~ column_label, scales = "free_y") +
    scale_colour_manual(values = pattern_colours) +
    scale_shape_manual(values = pattern_shapes) +
    labs(x = NULL, y = ylab, colour = "Effect pattern", shape = "Effect pattern")
  if (log) p + blank_125(d, y, c("P_label", "column_label"), c("sampler_label", "pattern_label")) +
    scale_y_log10(breaks = breaks_125, labels = labels_125) else p
}

save_plot(condition_plot("ess_per_second",
                         "Median cross-lag ESS per second", log = TRUE),
          "clpm_ess_per_second_by_sampler")
# The manuscript figure, at print size.
save_figure(condition_plot("ess_per_second", "Median ESS of the cross-lags per second", log = TRUE),
            "clpm_ess_per_second_by_sampler", height_cm = 15)
save_plot(condition_plot("run_time", "Warmup + sampling time (seconds)", log = TRUE),
          "clpm_runtime_by_sampler")
save_plot(condition_plot("rmse", "RMSE over cross-lagged effects"),
          "clpm_accuracy_rmse")

# PIP against posterior mean for every cross-lag, coloured by selection outcome.
d <- label_fits(cross)
d$outcome <- factor(ifelse(d$active, ifelse(d$pip >= 0.5, "True positive", "False negative"),
                           ifelse(d$pip >= 0.5, "False positive", "True negative")),
                    levels = c("True positive", "False negative", "False positive",
                               "True negative"))
p <- ggplot(d, aes(mean, pip, colour = outcome)) +
  geom_vline(xintercept = 0, colour = "grey80") +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey60") +
  geom_point(alpha = 0.6, size = 1.3) +
  facet_grid(P_label + N_label ~ sampler_label) +
  scale_colour_manual(values = c("True positive" = "#0072B2", "False negative" = "#CC79A7",
                                 "False positive" = "#D55E00", "True negative" = "grey60")) +
  labs(x = "Posterior mean", y = "Posterior inclusion probability", colour = NULL)
save_plot(p, "clpm_inclusion_probability_vs_mean", width = 10, height = 6)

# PIPs of each PDMP sampler against the NIMBLE reference.
d <- label_fits(agreement)
p <- ggplot(d, aes(ref_pip, pip)) +
  geom_abline(colour = "grey60") +
  geom_point(alpha = 0.4, size = 1, colour = "#0072B2") +
  facet_grid(prior_label ~ sampler_label) +
  coord_equal() +
  labs(x = "PIP, NIMBLE reference", y = "PIP, PDMP sampler")
save_plot(p, "clpm_pip_agreement_with_nimble", width = 8, height = 5.5)

message("Summaries: ", summary_dir, "\nFigures: ", plot_dir)

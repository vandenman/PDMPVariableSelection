# Appendix tables of the empirical RI-CLPM example (scripts/empirical_clpm.R).
#
# Run from the project root after the analysis:
#     Rscript scripts/table_empirical_clpm.R
# Writes two LaTeX tables for the Boomerang fits under the dependent prior (both
# chains averaged): the posterior inclusion probabilities of all cross-lags, and
# the posterior means of all lagged effects (autoregressive effects on the
# diagonal). Rows are predictors at wave t - 1, columns outcomes at wave t.

source("scripts/simulation_common.R")   # fit_parameters()

out_dir <- "results/empirical_liss_health"
table_dir <- "../proposal/tables"
dir.create(table_dir, showWarnings = FALSE)
labels <- c(anxious = "Anxious", down = "Down", calm_r = "Not calm",
            depressed = "Depressed", happy_r = "Not happy", poor_health = "Poor health",
            hinder_daily = "Hinders daily", hinder_social = "Hinders social",
            pain = "Pain", fatigue = "Fatigue", sleep = "Sleep problems")
variables <- names(labels)

effects <- do.call(rbind, lapply(1:2, function(chain) {
  f <- readRDS(file.path(out_dir, "fits", sprintf("dependent_AdaptiveBoomerang_chain%d.rds", chain)))
  p <- fit_parameters(f)
  ij <- regmatches(p$parameter, regexec("^A\\.([0-9]+)\\.([0-9]+)$", p$parameter))
  keep <- lengths(ij) == 3
  data.frame(outcome = variables[as.integer(sapply(ij[keep], `[`, 2))],
             predictor = variables[as.integer(sapply(ij[keep], `[`, 3))],
             pip = p$pip[keep], mean = p$mean[keep])
}))
effects <- aggregate(cbind(pip, mean) ~ predictor + outcome, effects, mean,
                     na.action = na.pass)    # the autoregressive effects have no PIP
effects$pip[effects$predictor == effects$outcome] <- NA   # always included

as_matrix <- function(value) {
  m <- matrix(NA_real_, length(variables), length(variables),
              dimnames = list(variables, variables))
  m[cbind(effects$predictor, effects$outcome)] <- effects[[value]]
  m
}

# A rotated-header LaTeX table; `bold` marks the cells to print in bold.
write_table <- function(m, digits, bold, caption, label, file) {
  cell <- ifelse(is.na(m), "", formatC(m, format = "f", digits = digits))
  cell <- sub("^(-?)0\\.", "\\1.", cell)                     # .12 instead of 0.12
  cell <- sub("^-(\\.0+)$", "\\1", cell)                       # no negative zeros
  cell <- sub("^-", "$-$", cell)                                # a minus sign, not a hyphen
  cell[bold & !is.na(m)] <- paste0("\\textbf{", cell[bold & !is.na(m)], "}")
  diag(cell)[is.na(diag(m))] <- "--"
  header <- paste0("\\rotatebox{90}{", labels[variables], "}", collapse = " & ")
  rows <- paste0(labels[variables], " & ", apply(cell, 1, paste, collapse = " & "), " \\\\")
  writeLines(c(
    "\\begin{table}[H]", "\\centering",
    paste0("\\caption{", caption, "}"), paste0("\\label{", label, "}"),
    "\\par\\medskip", "\\scriptsize", "\\setlength{\\tabcolsep}{3pt}",
    paste0("\\begin{tabular}{l", strrep("r", length(variables)), "}"),
    "\\toprule", paste0("Predictor & ", header, " \\\\"), "\\midrule",
    rows, "\\bottomrule", "\\end{tabular}", "\\end{table}"),
    file.path(table_dir, file))
}

pip <- as_matrix("pip")
write_table(pip, 2, pip > 0.5,
  paste("Posterior inclusion probabilities of the cross-lagged effects in the empirical example",
        "(dependent prior, Boomerang sampler, two chains).",
        "Rows are predictors and columns outcomes a year later.",
        "Inclusion probabilities above 0.50 are in bold."),
  "tab:empirical-pip", "empirical_pip.tex")
means <- as_matrix("mean")
write_table(means, 3, row(means) == col(means) | pip > 0.5,
  paste("Posterior means of the lagged effects in the empirical example",
        "(dependent prior, Boomerang sampler, two chains).",
        "Rows are predictors and columns outcomes a year later.",
        "The diagonal contains the autoregressive effects, which are always included;",
        "they and the cross-lagged effects with inclusion probability above 0.50 are in bold.",
        "The means of the cross-lagged effects are marginal means, averaged over the time an effect is excluded."),
  "tab:empirical-means", "empirical_means.tex")
cat(sprintf("%d cross-lags with PIP > 0.5; tables in %s\n", sum(pip > 0.5, na.rm = TRUE), table_dir))

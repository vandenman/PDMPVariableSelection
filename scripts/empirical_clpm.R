# Empirical example: a sparse RI-CLPM of mood and physical health in the LISS
# panel (11 variables, six annual waves 2008-2013, 2,670 persons).
#
# Run from the project root, after preparing the data with
# ../dataset/prepare_liss_health.R:
#     Rscript scripts/empirical_clpm.R
# The sections are: 1. settings  2. data  3. fits  4. summary  5. trajectory
# for the figure (scripts/plot_empirical_clpm.R).
#
# The selectable parameters are the 110 cross-lagged effects. Both priors of the
# simulation study are used (independent Bernoulli(0.2) inclusion; dependent
# Beta(1, 4) model prior with an exchangeable slab of correlation 0.25), with a
# slab sd of 0.1 instead of 0.35: within-person cross-lagged effects are small
# (benchmarks of .03, .07 and .12 for small, medium and large effects; Orth et
# al., 2024), and 0.1 keeps most prior mass on effects up to "large". Each
# prior is fitted with the Boomerang and Bouncy Particle samplers and with
# NIMBLE reversible jump, two chains each.

# The study script provides the model and fitting functions (fit_clpm()); its
# own simulation section is skipped.
worker <- identical(Sys.getenv("SIMULATION_WORKER"), "TRUE")
Sys.setenv(SIMULATION_WORKER = "TRUE")
source("scripts/simulation_clpm.R")
if (!worker) Sys.unsetenv("SIMULATION_WORKER")

# ---- 1. Settings --------------------------------------------------------------

panel_file <- "../dataset/liss_health_panel.csv"
slab_sd <- 0.1                 # slab_correlation (0.25) as in the simulation study
retained_T <- 40000            # retained PDMP trajectory (physical time); warmup
                               # is 250 P as in the simulation study
nimble_burnin <- 50000
nimble_iterations <- 50000
chains <- 2
workers <- as.integer(Sys.getenv("SIM_WORKERS", 12))
base_seed <- 20261001
out_dir <- "results/empirical_liss_health"

# ---- 2. Data ------------------------------------------------------------------

# One row per person, the waves side by side (all variables at wave 1, then at
# wave 2, ...). Each variable is standardised with its mean and sd over all
# person-waves, so that the cross-lags are on a common scale.
panel <- read.csv(panel_file)
variables <- setdiff(names(panel), c("id", "wave", "year"))
panel[variables] <- lapply(panel[variables], function(x) (x - mean(x)) / sd(x))
panel <- panel[order(panel$id, panel$wave), ]
n_waves <- length(unique(panel$wave))
z <- do.call(cbind, lapply(seq_len(n_waves), function(w)
  as.matrix(panel[panel$wave == w, variables])))
stopifnot(nrow(z) == length(unique(panel$id)), !anyNA(z))

# ---- 3. Fits ------------------------------------------------------------------

tasks <- expand.grid(chain = seq_len(chains), prior = c("dependent", "independent"),
                     sampler = c("AdaptiveBoomerang", "PreconditionedBPS", "nimble"),
                     stringsAsFactors = FALSE)
tasks$seed <- base_seed + seq_len(nrow(tasks))
tasks$id <- with(tasks, sprintf("%s_%s_chain%d", prior, sampler, chain))
# NIMBLE chains are the longest: start them first.
tasks <- tasks[order(tasks$sampler != "nimble"), ]

# A separate, shorter Boomerang run under the dependent prior that keeps its
# retained trajectory of the cross-lags (time step 0.05), for the figure.
trajectory_run <- function() {
  fit <- fit_clpm(z, nrow(z), length(variables), n_waves, "dependent",
                  "AdaptiveBoomerang", base_seed + 100, "trajectory",
                  retained = 500, trajectory_dt = 0.05)
  saveRDS(list(trajectory = fit$trajectory, variables = variables),
          file.path(out_dir, "trajectory.rds"))
}

fit_task <- function(task)
  fit_clpm(z, nrow(z), length(variables), n_waves, task$prior, task$sampler,
           task$seed, task$id)

if (!worker) {
  fits_dir <- run_tasks(tasks, setup_file = "scripts/empirical_clpm.R",
                        out_dir = out_dir, workers = workers)

  # ---- 4. Summary --------------------------------------------------------------

  results <- collect_fits(fits_dir)
  write.csv(results$fits, file.path(out_dir, "empirical_fits.csv"), row.names = FALSE)
  # Cross-lag A[i, j]: effect of variable j at wave t on variable i at wave t + 1
  # (PDMP output names it A.i.j, NIMBLE A[i, j]).
  p <- merge(results$parameters, results$fits[c("id", "prior", "sampler", "chain")])
  ij <- regmatches(p$parameter, regexec("^A\\W+([0-9]+)\\W+([0-9]+)\\W*$", p$parameter))
  p$outcome <- vapply(ij, function(m) if (length(m) == 3) variables[as.integer(m[2])] else NA_character_, "")
  p$predictor <- vapply(ij, function(m) if (length(m) == 3) variables[as.integer(m[3])] else NA_character_, "")
  cross <- p[!is.na(p$outcome) & p$outcome != p$predictor, ]
  write.csv(cross, file.path(out_dir, "empirical_cross_lags.csv"), row.names = FALSE)
  summary <- aggregate(cbind(pip, mean) ~ predictor + outcome + prior + sampler, cross, mean)
  summary <- summary[order(summary$prior, summary$sampler, -summary$pip), ]
  print(summary[summary$pip > 0.25, ], digits = 2, row.names = FALSE)

  # ---- 5. Trajectory for the figure --------------------------------------------

  trajectory_run()
}


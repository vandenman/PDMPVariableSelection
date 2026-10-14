# Simulation study for the random-intercept cross-lagged panel model (RI-CLPM).
#
# Run from the project root, either as a whole
#     Rscript scripts/simulation_clpm.R
# or step by step in an interactive session. The sections are:
#   1. settings          2. design          3. data generation
#   4. fitting one task  5. running all tasks and collecting the results
#
# The selectable parameters are the P(P - 1) cross-lagged effects, for P = 4, 8
# and 12 variables (12, 56 and 132 cross-lags). Two priors
# are compared: an independent prior (inclusion probability 0.2, Gaussian slab
# with sd 0.35) and a dependent prior (Beta(1, 4) model prior, exchangeable
# Gaussian slab with sd 0.35 and correlation 0.25). The PDMP samplers use the
# sufficient-statistic likelihood; NIMBLE provides a reference for both priors,
# with the same (sufficient-statistic) likelihood and reversible jump.

source("scripts/simulation_common.R")

# ---- 1. Settings --------------------------------------------------------------

Ps <- c(4, 8, 12)              # numbers of variables
reps <- 10                     # replications per data condition
reps_first <- 5                # replications 6-10 of P = 12 run last (see below)
workers <- as.integer(Sys.getenv("SIM_WORKERS", 10))  # parallel R processes
warmup_per_variable <- 250     # PDMP warmup is 250 * P (physical time)
retained_T <- 4000             # retained PDMP trajectory (physical time)
nimble_burnin <- 2000
nimble_iterations <- 4000      # retained NIMBLE iterations
slab_sd <- 0.1                  # slab sd of the cross-lags, as in the empirical example
slab_correlation <- 0.25
base_seed <- 20260702
out_dir <- "simulation_results/clpm_study"

# ---- 2. Design ----------------------------------------------------------------

# The two sample sizes bracket the empirical example (N = 2,670): cross-lags of
# 0.07 are hard to detect with 1,000 persons and clear with 5,000.
conditions <- expand.grid(N = c(1000, 5000), P = Ps, T = 4,
                          pattern = c("null", "sender", "mixed"),
                          rep = seq_len(reps), stringsAsFactors = FALSE)
conditions$data_seed <- base_seed + seq_len(nrow(conditions))

methods <- expand.grid(
  prior = c("independent", "dependent"),
  sampler = c("PreconditionedZigZag", "PreconditionedBPS", "AdaptiveBoomerang",
              "nimble"), stringsAsFactors = FALSE)
methods$method_no <- seq_len(nrow(methods))

tasks <- merge(conditions, methods, by = NULL)
tasks$id <- with(tasks, sprintf("N%d_P%d_%s_rep%02d__%s_%s", N, P, pattern, rep,
                                prior, sampler))
# Seeds depend only on the condition and method, so raising `reps` later keeps
# the seeds (and saved results) of the existing replications.
tasks$sampler_seed <- tasks$data_seed * 100 + tasks$method_no
# Longest fits first, so that they do not end up in the tail of the run: larger
# P first, and within P Zig-Zag, NIMBLE, Bouncy Particle, Boomerang. The P = 12
# replications after the first `reps_first` run at the very end, so a run that
# is cut short still has the full design with `reps_first` replications at P = 12.
# Ties are run in random (but reproducible) order.
set.seed(base_seed)
tasks <- tasks[order(tasks$P == 12 & tasks$rep > reps_first, -tasks$P,
                     match(tasks$sampler, c("PreconditionedZigZag", "nimble",
                                            "PreconditionedBPS", "AdaptiveBoomerang")),
                     runif(nrow(tasks))), ]
rownames(tasks) <- NULL

# ---- 3. Data generation -------------------------------------------------------

# Cross-lagged matrix A (row = outcome, column = predictor): autoregressive
# effects 0.35 on the diagonal and cross-lags of size `signal` (0.07, the
# medium within-person effect of Orth et al., 2024) in the pattern
#   null    no cross-lags
#   sender  variable 1 predicts all other variables
#   mixed   a random quarter of the cross-lags, with alternating signs
clpm_transition <- function(P, pattern, signal = 0.07) {
  A <- diag(0.35, P)
  if (pattern == "sender") A[-1, 1] <- signal
  if (pattern == "mixed") {
    cross <- which(row(A) != col(A))
    active <- sample(cross, max(2, round(0.25 * length(cross))))
    A[active] <- signal * sample(rep(c(-1, 1), length.out = length(active)))
  }
  A
}

# RI-CLPM data: y[t] = eta + w[t], with w[t] = A w[t - 1] + innovation.
simulate_clpm <- function(N, P, T, pattern, seed) {
  set.seed(seed)
  A <- clpm_transition(P, pattern)
  compound <- function(variance, covariance) {
    S <- matrix(covariance, P, P); diag(S) <- variance; S
  }
  draw <- function(S) matrix(stats::rnorm(N * P), N) %*% chol(S)
  eta <- draw(compound(0.60^2, 0.08))
  w <- draw(compound(0.70^2, 0.08))
  z <- eta + w
  for (t in 2:T) {
    w <- w %*% t(A) + draw(compound(0.50^2, 0.04))
    z <- cbind(z, eta + w)
  }
  list(z = z, A = A)
}

# Stan data for the sufficient-statistic model, with the slab of `prior`.
clpm_stan_data <- function(z, N, P, T, prior) {
  E <- P * (P - 1)
  rho <- if (prior == "dependent") slab_correlation else 0
  list(N = N, P = P, T = T, prior_only = 0L, z_sum = colSums(z),
       z_crossprod = crossprod(z), slab_mean = rep(0, E),
       slab_cov = slab_sd^2 * ((1 - rho) * diag(E) + rho))
}

# ---- 4. Fitting one task ------------------------------------------------------

# One fit of the RI-CLPM to data z (N persons, columns are wave blocks of the P
# variables) with the given prior and sampler. Also used by the empirical
# example (scripts/empirical_clpm.R), which also uses `retained` and
# `trajectory_dt` (see fit_pdmp()) for a figure of the trajectory.
fit_clpm <- function(z, N, P, T, prior, sampler, seed, id, retained = retained_T,
                     trajectory_dt = NULL) {
  task <- list(N = N, P = P, T = T, prior = prior, sampler = sampler,
               sampler_seed = seed, id = id)
  data <- list(z = z)
  stan_data <- clpm_stan_data(data$z, task$N, task$P, task$T, task$prior)
  # Every method starts at the posterior mode with all cross-lags included.
  mode <- stan_mode("clpm_optimized", stan_data)
  if (task$sampler == "nimble") {
    if (!exists("run_clpm_nimble_rj")) source("scripts/simulation_nimble.R")
    fit <- run_clpm_nimble_rj(
      list(z = data$z, P = task$P, T = task$T),
      niter = nimble_burnin + nimble_iterations, nburnin = nimble_burnin,
      prior_type = task$prior, slab_sd = slab_sd,
      slab_correlation = slab_correlation, seed = task$sampler_seed,
      # New cross-lags are proposed from, and moved on, the scale of the slab.
      proposal_scale = slab_sd,
      inits = clpm_nimble_inits_from_mode(mode))
    s <- fit$fit$posterior_summary
    result <- list(
      parameters = data.frame(parameter = s$parameter, selectable = TRUE,
        mean = s$posterior_mean, sd = NA_real_, pip = s$inclusion_probability,
        ess_mean = s$ess, ess_pip = NA_real_, transitions = NA_real_),
      timing = c(total = fit$fit$stats$elapsed_time + fit$fit$stats$compile_time,
                 warmup = NA_real_, sampling = fit$fit$stats$elapsed_time),
      # Retained draws of the selectable effects, for later recomputation.
      samples = fit$fit$target_samples,
      convergence = nimble_timed_batches(fit$fit$batches,
        colnames(fit$fit$target_samples), nimble_burnin,
        fit$fit$stats$elapsed_time))
  } else {
    prior_standata <- stan_data
    prior_standata$prior_only <- 1L
    parameter_names <- stan_names("clpm_optimized", stan_data)
    ij <- regmatches(parameter_names,
                     regexec("^A\\.([0-9]+)\\.([0-9]+)$", parameter_names))
    selectable <- vapply(ij, function(m) length(m) == 3 && m[2] != m[3], TRUE)
    sticky_prior <- if (task$prior == "independent") {
      list(model_prior = bernoulli(0.2),
           parameter_prior = rep(stats::dnorm(0, sd = slab_sd),
                                 length(parameter_names)))
    } else {
      list(model_prior = betabernoulli(1, 4),
           slab_prior = exchangeable_gaussian_slab(
             mean = 0, u = slab_sd^2 * (1 - slab_correlation),
             v = slab_sd^2 * slab_correlation,
             coef = parameter_names[selectable]))
    }
    # The initial speed of the cross-lags comes from the Laplace standard
    # deviations at the mode, as in the OMRF study.
    result <- fit_pdmp(
      "clpm_optimized", stan_data, prior_standata, task$sampler, sticky_prior,
      selectable = selectable, x0 = mode$mode, x0_sd = mode$sd,
      warmup_T = warmup_per_variable * task$P,
      retained_T = retained, seed = task$sampler_seed, trajectory_dt = trajectory_dt,
      work_dir = file.path(tempdir(), task$id))
    unlink(file.path(tempdir(), task$id), recursive = TRUE)
  }
  result
}

fit_task <- function(task) {
  data <- simulate_clpm(task$N, task$P, task$T, task$pattern, task$data_seed)
  result <- fit_clpm(data$z, task$N, task$P, task$T, task$prior, task$sampler,
                     task$sampler_seed, task$id)
  # True cross-lag A[i, j] for parameter "A.i.j" (NIMBLE writes "A[i, j]").
  ij <- regmatches(result$parameters$parameter, regexec(
    "^A\\W+([0-9]+)\\W+([0-9]+)\\W*$", result$parameters$parameter))
  result$parameters$true_value <- vapply(ij, function(m)
    if (length(m) == 3) data$A[as.integer(m[2]), as.integer(m[3])] else NA_real_,
    numeric(1))
  result
}

# ---- 5. Run all tasks and collect the results ---------------------------------

# The parallel workers source this file for its definitions only; run_tasks()
# marks them with SIMULATION_WORKER so that they skip this section.
if (!identical(Sys.getenv("SIMULATION_WORKER"), "TRUE")) {
  fits_dir <- run_tasks(tasks, setup_file = "scripts/simulation_clpm.R",
                        out_dir = out_dir, workers = workers)
  results <- collect_fits(fits_dir)
  write.csv(results$fits, file.path(out_dir, "clpm_fits.csv"), row.names = FALSE)
  write.csv(results$parameters, file.path(out_dir, "clpm_parameters.csv"),
            row.names = FALSE)
}

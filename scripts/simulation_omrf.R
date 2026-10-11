# Simulation study for the ordinal Markov random field (OMRF).
#
# Run from the project root, either as a whole
#     Rscript scripts/simulation_omrf.R
# or step by step in an interactive session. The sections are:
#   1. settings          2. design          3. data generation
#   4. fitting one task  5. running all tasks and collecting the results
#
# Two comparisons share the data sets:
#   gradient  independent prior (Beta(1, 4) model prior, standard Cauchy slab),
#             full versus subsampled gradients for Zig-Zag, Bouncy Particle and
#             Boomerang;
#   prior     independent versus node-specific Gaussian scale slab, full and
#             subsampled gradients for Zig-Zag, Bouncy Particle and Boomerang.
# NIMBLE provides a reference for both priors in every condition.

source("scripts/simulation_common.R")

# ---- 1. Settings --------------------------------------------------------------

reps <- 10                     # replications per data condition
workers <- as.integer(Sys.getenv("SIM_WORKERS", 30))  # parallel R processes
warmup_per_node <- 20          # PDMP warmup is 20 * P (physical time), started at the
                               # posterior mode so the model size settles quickly
retained_T <- 1000             # retained PDMP trajectory (physical time)
nimble_burnin <- 2000
nimble_iterations <- 8000      # retained NIMBLE iterations
slab_base_scale <- 1.25        # node prior: slab sd = 1.25 sqrt(scale_i scale_j), so
                               # that its density at zero equals the Cauchy(0, 1) slab
node_scale_sd <- 0.5           # node prior: node scales ~ lognormal(0, 0.5)
base_seed <- 20260702
out_dir <- "simulation_results/omrf_study"

# ---- 2. Design ----------------------------------------------------------------

conditions <- expand.grid(N = c(200, 500), P = c(10, 20, 30), K = 4,
                          structure = c("random", "hub"), rep = seq_len(reps),
                          stringsAsFactors = FALSE)
conditions$data_seed <- base_seed + seq_len(nrow(conditions))

dynamics <- c("PreconditionedZigZag", "PreconditionedBPS", "AdaptiveBoomerang")
methods <- rbind(
  data.frame(comparison = "gradient", prior = "independent",
             sampler = rep(dynamics, 2),
             gradient = rep(c("full", "subsampled"), each = 3)),
  data.frame(comparison = "prior", prior = "node",
             sampler = c(dynamics, "AdaptiveBoomerang"),
             gradient = c("full", "full", "full", "subsampled")),
  data.frame(comparison = "reference", prior = c("independent", "node"),
             sampler = "nimble", gradient = NA),
  # Added last, so that the method numbers (and seeds) of the others are kept.
  data.frame(comparison = "prior", prior = "node",
             sampler = c("PreconditionedBPS", "PreconditionedZigZag"),
             gradient = "subsampled"))
methods$method_no <- seq_len(nrow(methods))

tasks <- merge(conditions, methods, by = NULL)
tasks$id <- with(tasks, sprintf("N%d_P%d_%s_rep%02d__%s_%s_%s", N, P, structure,
                                rep, prior, sampler, gradient))
# Seeds depend only on the condition and method, so raising `reps` later keeps
# the seeds (and saved results) of the existing replications.
tasks$sampler_seed <- tasks$data_seed * 100 + tasks$method_no
# Longest fits first, so that they do not end up in the tail of the run: larger
# P first, and within P the node prior and the Zig-Zag sampler (most events).
# Ties are run in random (but reproducible) order.
set.seed(base_seed)
tasks <- tasks[order(-tasks$P, tasks$prior != "node",
                     tasks$sampler != "PreconditionedZigZag",
                     runif(nrow(tasks))), ]
rownames(tasks) <- NULL

# ---- 3. Data generation -------------------------------------------------------

# True network: a fraction `density` of the P(P - 1)/2 edges is active. "random"
# chooses the active edges uniformly; "hub" picks ceiling(P / 5) random hub
# nodes and chooses the active edges among the edges that touch a hub. Edge
# magnitudes vary by node: |beta_ij| = signal * sqrt(a_i a_j) with node
# strengths a_i = exp(0.5 z_i), z_i ~ N(0, 1); signs are random.
# Thresholds give every node the category marginals `marginals` in the absence
# of interactions (alpha_u = log p_u / p_0), and are centred on the mean field:
# subtracting u * sum_j beta_ij (K - 1) / 2 makes the logits
# alpha_u + u * sum_j beta_ij (x_j - (K - 1) / 2), so that nodes with several
# positive edges are not pushed into one extreme category.
omrf_truth <- function(P, K, structure, seed, density = 0.2, signal = 0.25,
                       strength_sd = 0.5, marginals = c(0.15, 0.35, 0.35, 0.15)) {
  set.seed(seed)
  edges <- edge_index(P)
  n_active <- round(density * nrow(edges))
  candidates <- edges$edge
  hubs <- integer(0)
  if (structure == "hub") {
    hubs <- sample(P, ceiling(P / 5))
    candidates <- edges$edge[edges$i %in% hubs | edges$j %in% hubs]
  }
  active <- sort(candidates[sample.int(length(candidates), n_active)])
  signs <- sample(c(-1, 1), n_active, replace = TRUE)
  strength <- exp(strength_sd * stats::rnorm(P))
  i <- edges$i[active]
  j <- edges$j[active]
  interactions <- matrix(0, P, P)
  interactions[cbind(i, j)] <- signs * signal * sqrt(strength[i] * strength[j])
  interactions <- interactions + t(interactions)
  alpha <- log(marginals[-1] / marginals[1])
  thresholds <- t(replicate(P, alpha + stats::rnorm(K - 1, sd = 0.1))) -
    outer(rowSums(interactions) * (K - 1) / 2, seq_len(K - 1))
  list(interactions = interactions, thresholds = thresholds,
       edge_values = interactions[cbind(edges$i, edges$j)],
       strength = strength, hubs = hubs)
}

# N responses from the OMRF, simulated with bgms (1000 Gibbs iterations). bgms
# uses the same parameterization: pairwise = interactions, main = thresholds.
simulate_omrf <- function(N, P, K, structure, seed) {
  truth <- omrf_truth(P, K, structure, seed)
  X <- bgms::simulate_mrf(num_states = N, num_variables = P,
                          num_categories = K - 1, pairwise = truth$interactions,
                          main = truth$thresholds, iter = 1000, seed = seed)
  storage.mode(X) <- "integer"
  stan_data <- list(N = N, P = P, K = K, prior_only = 0L, X = X,
                    seen = rep(K, P), prior_cauchy_scale = 1,
                    prior_threshold_alpha = 1, prior_threshold_beta = 1,
                    prior_node_scale_sd = node_scale_sd,
                    slab_base_scale = slab_base_scale)
  list(stan_data = stan_data, truth = truth)
}

# ---- 4. Fitting one task ------------------------------------------------------

# Sticky prior for the PDMP sampler: Beta(1, 4) model prior with either the
# independent standard Cauchy slab (in the Stan model) or the node-specific
# Gaussian scale slab, whose log sd is log(slab_base_scale) + (log node scale i
# + log node scale j) / 2.
omrf_sticky_prior <- function(prior, P, n_parameters) {
  if (prior == "independent") {
    return(list(model_prior = betabernoulli(1, 4),
                parameter_prior = rep(stats::dcauchy(0), n_parameters)))
  }
  edges <- edge_index(P)
  design <- Matrix::sparseMatrix(
    i = rep(edges$edge, each = 2), j = c(rbind(edges$i, edges$j)),
    x = 0.5, dims = c(nrow(edges), P))
  list(model_prior = betabernoulli(1, 4),
       slab_prior = loglinear_gaussian_scale_slab(
         log_base_sd = log(slab_base_scale), logscale = "node_scale",
         logscale_design = design, coef = "interactions_0"))
}

# Subsampled gradient. The gradient is estimated from a random group of persons
# with a Hessian control variate around the nearest of eight anchors, which
# fit_pdmp() takes from full-gradient sampling after the warmup; "nearest" is
# measured in the metric of the mean anchor Hessian. Bouncy
# Particle with the independent prior uses fixed groups of 16 persons with joint
# signed event-rate bounds; the other dynamics and priors use one person per
# estimate.
omrf_subsampling <- function(stan_data, prior_standata, n_parameters, sampler,
                             prior) {
  grouped <- sampler == "PreconditionedBPS" && prior == "independent"
  # The envelope of the paper (per-person direct-factor roof with Hessian
  # control variates, analytic prior) is the only one PDMPSamplersR provides.
  # Versions before its October 2026 clean-up offered several, with other
  # defaults, and would silently build a different envelope from this call.
  if ("geometry" %in% names(formals(omrf_residual_envelope)))
    stop("this script needs PDMPSamplersR after its October 2026 clean-up")
  envelope <- omrf_residual_envelope(
    X = stan_data$X, seen = stan_data$seen, thresholds = "thresholds_0",
    interactions = "interactions_0", exact_selected_refinement = TRUE,
    person_signed_roof = grouped, signed_only_person_roof = grouped,
    joint_signed_group_size = if (grouped) 16L else 0L,
    shared_current_categorical = grouped)
  list(subsampling = stan_subsampling(
         n_observations = stan_data$N, subsample_size = if (grouped) 16L else 1L,
         prior_standata = prior_standata, residual_envelope = envelope,
         anchor = rep(0, n_parameters)),
       args = list(subsampling_use_anchor_bank = TRUE,
                   subsampling_anchor_bank_capacity = 8L,
                   subsampling_require_no_main_full_data = TRUE),
       bank_metric = "hessian")
}

fit_task <- function(task) {
  data <- simulate_omrf(task$N, task$P, task$K, task$structure, task$data_seed)
  stan_data <- data$stan_data
  # Every method starts at the posterior mode of the independent-slab model
  # with every edge included; under the node prior the node scales start at 1.
  # With task$start_distance (scripts/simulation_far_start.R), every method
  # starts that far from the mode, and the PDMP samplers get no information
  # from the mode (unit scales, Boomerang reference centred at zero).
  mode <- stan_mode("omrf", stan_data)
  start_distance <- task$start_distance %||% 0
  if (start_distance > 0) mode <- far_start(mode, start_distance, task$data_seed)
  if (task$sampler == "nimble") {
    if (!exists("run_omrf_nimble_rj")) source("scripts/simulation_nimble.R")
    prior_type <- paste0("beta_binomial_", task$prior)
    fit <- run_omrf_nimble_rj(stan_data, prior_type = prior_type,
      burnin = nimble_burnin, iterations = nimble_iterations,
      seed = task$sampler_seed,
      inits = omrf_nimble_inits_from_mode(mode, stan_data, prior_type))
    s <- fit$fit$posterior_summary
    result <- list(
      parameters = data.frame(parameter = s$parameter, selectable = TRUE,
        mean = s$posterior_mean, sd = NA_real_, pip = s$inclusion_probability,
        ess_mean = s$ess, ess_pip = NA_real_, transitions = NA_real_),
      timing = c(total = fit$fit$stats$workflow_time, warmup = NA_real_,
                 sampling = fit$fit$stats$elapsed_time),
      # Retained draws of the edges (zero when excluded) and, for the node
      # prior, of the log scales, for later recomputation.
      samples = list(interactions = fit$fit$samples,
                     log_scales = fit$fit$scale_samples),
      convergence = nimble_timed_batches(fit$fit$batches,
        paste0("interactions[", seq_len(ncol(fit$fit$samples)), "]"),
        nimble_burnin, fit$fit$stats$elapsed_time, colnames(fit$fit$samples)))
  } else {
    model <- if (task$prior == "node") "omrf_dependent" else "omrf"
    # Prior-only data; the fields present tell the sampler which prior is used.
    prior_standata <- stan_data
    prior_standata$prior_only <- 1L
    if (task$prior == "node") prior_standata$prior_cauchy_scale <- NULL else
      prior_standata[c("prior_node_scale_sd", "slab_base_scale")] <- NULL
    parameter_names <- stan_names(model, stan_data)
    selectable <- grepl("^interactions_0\\.", parameter_names)
    subsampled <- identical(task$gradient, "subsampled")
    x0 <- numeric(length(parameter_names))
    x0_sd <- rep(1, length(parameter_names))
    x0[seq_along(mode$mode)] <- mode$mode  # thresholds and interactions first
    x0_sd[seq_along(mode$sd)] <- mode$sd
    if (start_distance > 0) x0_sd <- NULL
    result <- fit_pdmp(
      model, stan_data, prior_standata, task$sampler,
      sticky_prior = omrf_sticky_prior(task$prior, task$P,
                                       length(parameter_names)),
      selectable = selectable, x0 = x0, x0_sd = x0_sd,
      warmup_T = warmup_per_node * task$P, retained_T = retained_T,
      seed = task$sampler_seed,
      work_dir = file.path(tempdir(), task$id),
      subsampling = if (subsampled) function(names) omrf_subsampling(
        stan_data, prior_standata, length(names), task$sampler, task$prior))
    unlink(file.path(tempdir(), task$id), recursive = TRUE)
  }
  edges <- edge_index(task$P)
  result$parameters$true_value <- data$truth$edge_values[match(
    result$parameters$parameter, paste0("interactions_0.", edges$edge))]
  result
}

# ---- 5. Run all tasks and collect the results ---------------------------------

# The parallel workers source this file for its definitions only; run_tasks()
# marks them with SIMULATION_WORKER so that they skip this section.
if (!identical(Sys.getenv("SIMULATION_WORKER"), "TRUE")) {
  fits_dir <- run_tasks(tasks, setup_file = "scripts/simulation_omrf.R",
                        out_dir = out_dir, workers = workers)
  results <- collect_fits(fits_dir)
  write.csv(results$fits, file.path(out_dir, "omrf_fits.csv"), row.names = FALSE)
  write.csv(results$parameters, file.path(out_dir, "omrf_parameters.csv"),
            row.names = FALSE)
}

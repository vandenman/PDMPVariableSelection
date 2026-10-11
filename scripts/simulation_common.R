# Shared code for the simulation studies (simulation_omrf.R, simulation_clpm.R).
#
# Main functions:
#   fit_pdmp()          one sticky PDMP fit, summarised from online statistics
#   online_estimates()  posterior summaries from those online statistics
#   run_tasks()         run a table of fits in parallel, one file per fit
#   collect_fits()      combine the saved fits into two tables
#
# Every PDMP fit uses streaming trace storage with online summaries: the
# sampler accumulates exact continuous-time batch integrals instead of storing
# the event trace, so memory use does not grow with the number of events.

library(PDMPSamplersR)
# During development, load PDMPSamplersR from a source checkout instead.
if (nzchar(Sys.getenv("PDMPSAMPLERSR_SOURCE")))
  pkgload::load_all(Sys.getenv("PDMPSAMPLERSR_SOURCE"), quiet = TRUE)

`%||%` <- function(x, y) if (is.null(x)) y else x

stan_file <- function(model) file.path("stanmodels", paste0(model, ".stan"))

# Names of the unconstrained parameters of a Stan model, in sampler order.
stan_names <- function(model, stan_data) {
  compiled <- compile_pdmp_stan_model(normalizePath(stan_file(model)))
  stan_parameter_mapping(compiled, stan_data)$name
}

# Posterior mode of a Stan model on the unconstrained scale, including the
# log-Jacobian (the mode cmdstanr's laplace() uses), found by L-BFGS on
# BridgeStan's log density, and the Laplace standard deviations there (from
# the inverse Hessian, or 1 where the Hessian is not positive definite).
# The mode is the starting point of the PDMP warmup; the standard deviations
# set the initial speed of the selectable coordinates (see fit_pdmp()).
stan_mode <- function(model, stan_data) {
  compiled <- compile_pdmp_stan_model(normalizePath(stan_file(model)))
  data_file <- tempfile(fileext = ".json")
  on.exit(unlink(data_file), add = TRUE)
  write_stan_json(stan_data, data_file)
  JuliaCall::julia_command("using BridgeStan")
  JuliaCall::julia_command(sprintf(
    'mode_model = BridgeStan.StanModel("%s", "%s"; warn = false);',
    compiled, data_file))
  JuliaCall::julia_command(paste(
    "mode_nlp(x) = -BridgeStan.log_density(mode_model, x; jacobian = true);",
    "mode_grad(x) = -BridgeStan.log_density_gradient(mode_model, x;",
    "jacobian = true)[2];",
    "mode_hess(x) = -BridgeStan.log_density_hessian(mode_model, x;",
    "jacobian = true)[3];"))
  d <- JuliaCall::julia_eval("Int(BridgeStan.param_unc_num(mode_model))")
  mode <- stats::optim(numeric(d),
    function(x) JuliaCall::julia_call("mode_nlp", x),
    function(x) JuliaCall::julia_call("mode_grad", x),
    method = "L-BFGS-B", control = list(maxit = 5000))$par
  H <- JuliaCall::julia_call("mode_hess", mode)
  sd <- tryCatch(sqrt(diag(chol2inv(chol((H + t(H)) / 2)))),
                 error = function(e) rep(1, d))
  # The mode on the constrained scale, as arrays named after the Stan
  # parameters, which gives the NIMBLE reference the same starting point.
  JuliaCall::julia_assign("mode_x", mode)
  values <- JuliaCall::julia_eval(
    "BridgeStan.param_constrain(mode_model, mode_x)")
  names(values) <- JuliaCall::julia_eval(
    "BridgeStan.param_names(mode_model)")
  list(mode = mode, sd = sd, constrained = stan_arrays(values))
}

# A start away from the posterior mode (scripts/simulation_far_start.R): every
# unconstrained coordinate is moved `distance` from the mode in a random
# direction (seeded by `seed`, the same for all methods), and the constrained
# values used by NIMBLE are updated. Call directly after stan_mode(), whose
# BridgeStan model it uses.
far_start <- function(mode, distance, seed) {
  set.seed(seed)
  x <- mode$mode + distance * sample(c(-1, 1), length(mode$mode), replace = TRUE)
  JuliaCall::julia_assign("mode_x", x)
  values <- JuliaCall::julia_eval("BridgeStan.param_constrain(mode_model, mode_x)")
  names(values) <- JuliaCall::julia_eval("BridgeStan.param_names(mode_model)")
  mode$mode <- x
  mode$constrained <- stan_arrays(values)
  mode
}

# Named vector with Stan's flattened names ("A.1.2") to a list of arrays.
stan_arrays <- function(values) {
  block <- sub("\\..*$", "", names(values))
  lapply(split(values, factor(block, unique(block))), function(v) {
    idx <- do.call(rbind, lapply(strsplit(names(v), ".", fixed = TRUE),
                                 function(s) as.integer(s[-1])))
    if (is.null(idx) || ncol(idx) == 0) return(unname(v))
    if (ncol(idx) == 1) return(unname(v[order(idx[, 1])]))
    out <- array(NA_real_, apply(idx, 2, max))
    out[idx] <- v
    out
  })
}

# All pairs i < j of p nodes, in the order of the Stan interaction vector.
edge_index <- function(p) {
  idx <- utils::combn(seq_len(p), 2)
  data.frame(edge = seq_len(ncol(idx)), i = idx[1, ], j = idx[2, ])
}

# Settings shared by every PDMP fit (see the supplement for their meaning).
pdmp_settings <- list(
  algorithm = "GridThinningStrategy",
  grid_n = 5L, grid_t_max = 4, grid_bound = "value_quadratic",
  grid_curvature_bound = 10000, curvature_backend = "exact",
  c0 = 0.01, adaptive_scheme = "diagonal", post_warmup_simplify = TRUE,
  lazy_low_tightness_threshold = 0, lazy_max_low_tightness_rejections = 100000L,
  lazy_max_rejections = 100000L,
  online_batches = 400L,         # fine batches; estimates use 40 of them
  warmup_batches = 100L,         # warmup batches, for the time to convergence
  anchor_bank_T = 20,            # full-gradient time that supplies the anchors
  selectable_speed = 3,          # initial speed of selectable coordinates, in
                                 # Laplace standard deviations (see fit_pdmp)
  grid_points = 100              # retained-phase state snapshots, for traces
)

# Velocity refreshment. Bouncy Particle needs enough refreshment to mix; the
# Boomerang default is kept; Zig-Zag has none.
refresh_rate <- c(PreconditionedZigZag = 0, PreconditionedBPS = 1,
                  AdaptiveBoomerang = 0.1)

# Warmup metric updates, every this much physical time (NA: every warmup_T / 10).
# Bouncy Particle mixes poorly in model space unless its metric is updated
# often; for Zig-Zag frequent updates make the warmup very expensive.
adaptation_interval <- c(PreconditionedZigZag = NA, PreconditionedBPS = 1,
                         AdaptiveBoomerang = NA)

# One sticky PDMP fit.
#   stan_data        data for the Stan model `model`
#   prior_standata   the same data with prior_only = 1
#   sticky_prior     list(model_prior, parameter_prior, slab_prior) for
#                    PDMPSamplersR's sticky sampler
#   selectable       logical, which unconstrained coordinates can stick at zero
#   x0               initial position (unconstrained scale)
#   x0_sd            optional Laplace standard deviations at x0 (stan_mode()).
#                    The flow starts from this Laplace approximation: the other
#                    coordinates get scale x0_sd and, for the Boomerang, reference
#                    centre x0, which warmup adaptation then refines. Selectable
#                    coordinates are centred at zero and are not adapted (a stuck
#                    coordinate carries no information about its scale), so they
#                    keep their initial speed of `selectable_speed` standard
#                    deviations; the default 1 is far too fast for small effects.
#                    Without x0_sd, every coordinate starts with scale 1 and,
#                    for the Boomerang, reference centre 0 (no information from
#                    the mode).
#   full_gradient    optional analytic full-gradient provider (OMRF)
#   subsampling      optional function(parameter_names) returning
#                    list(subsampling = stan_subsampling(...), args = extra
#                    pdmp_sample_from_stanmodel arguments, bank_metric =
#                    optional anchor metric of position_anchor_bank())
#   work_dir         directory for the online summary files
#   trajectory_dt    optional time step: record the state of the retained
#                    trajectory at this spacing (instead of 100 snapshots) and
#                    return the selectable coefficients (for figures)
#   curvature_backend "exact" (Stan Hessian-vector products) or
#                    "finite_difference" (shifted gradients) for the grid
#                    thinning curvature of full-gradient runs
# Returns list(parameters = data.frame, timing = named numeric, counts = ...).
# timing: `warmup` and `sampling` are the sampler's own phase timers, so they
# exclude Julia compilation and model loading; for subsampled fits `warmup` also
# includes the full-gradient anchor stage and building the anchor bank. `total`
# is the wall time of the sampler calls, compilation included.
#
# Without subsampling this is one run: adaptive warmup, then the retained
# trajectory. With subsampling it is two runs. The first is the same
# full-gradient warmup, followed by `anchor_bank_T` of full-gradient sampling
# with the adapted metric; its states give the frozen anchors of the gradient
# control variates. The second continues from its endpoint with the same metric
# and subsampled gradients for the retained trajectory.
fit_pdmp <- function(model, stan_data, prior_standata, dynamic, sticky_prior,
                     selectable, x0, warmup_T, retained_T, seed, work_dir,
                     x0_sd = NULL, full_gradient = NULL, subsampling = NULL,
                     curvature_backend = pdmp_settings$curvature_backend,
                     trajectory_dt = NULL) {
  compiled <- compile_pdmp_stan_model(normalizePath(stan_file(model)))
  parameter_names <- stan_parameter_mapping(compiled, stan_data)$name
  common <- c(pdmp_settings[c("algorithm", "grid_n", "grid_t_max", "grid_bound",
    "grid_curvature_bound", "curvature_backend", "c0", "adaptive_scheme",
    "post_warmup_simplify", "lazy_low_tightness_threshold",
    "lazy_max_low_tightness_rejections", "lazy_max_rejections")], list(
    path_to_stanmodel = compiled, standata = stan_data, flow = dynamic,
    refresh_rate = refresh_rate[[dynamic]],
    warmup_adaptation_interval = if (!is.na(adaptation_interval[[dynamic]]))
      adaptation_interval[[dynamic]],
    sticky = TRUE, can_stick = selectable,
    model_prior = sticky_prior$model_prior,
    parameter_prior = sticky_prior$parameter_prior,
    slab_prior = sticky_prior$slab_prior,
    prior_stanmodel = if (is.null(sticky_prior$slab_prior)) NULL else compiled,
    show_progress = FALSE, seed = seed, materialize = FALSE,
    trace_storage = "streaming"))
  common$curvature_backend <- curvature_backend
  # One sampler call whose retained phase is summarised online in `stage`.
  run <- function(stage, grid_spacing = NA, ...) {
    dir <- file.path(work_dir, stage)
    unlink(dir, recursive = TRUE)
    args <- c(common, list(...), list(streaming_trace_directory = dir,
      online_summaries = online_summary_options(pdmp_settings$online_batches,
        grid_spacing, warmup_batches = pdmp_settings$warmup_batches)))
    args <- args[!duplicated(names(args), fromLast = TRUE)]  # stage overrides
    elapsed <- system.time(fit <- do.call(pdmp_sample_from_stanmodel, args))
    list(fit = fit, online = read_online_summaries(dir),
         online_warmup = tryCatch(read_online_summaries(dir, "warmup"),
                                  error = function(e) NULL),
         elapsed = elapsed[["elapsed"]])
  }
  full_gradient_run <- function(stage, retained, grid_spacing = NA) run(
    stage, grid_spacing, full_gradient = full_gradient, t0 = 0,
    t_warmup = warmup_T, T = warmup_T + retained, x0 = x0,
    flow_cov = if (is.null(x0_sd)) diag(length(x0)) else diag(ifelse(selectable,
      (pdmp_settings$selectable_speed * x0_sd)^2, x0_sd^2)),
    flow_mean = if (is.null(x0_sd)) numeric(length(x0)) else ifelse(selectable, 0, x0),
    flow_specification = "initial",
    # The prior-only data are needed for a separate slab.
    prior_standata = if (!is.null(sticky_prior$slab_prior)) prior_standata,
    support_boundary = support_boundary_control(
      mode = if (dynamic == "PreconditionedBPS")
        "line_search_truncated_refresh" else "error"))

  main_grid <- if (!is.null(trajectory_dt)) trajectory_dt else
    if (pdmp_settings$grid_points > 0) retained_T / pdmp_settings$grid_points else NA
  if (is.null(subsampling)) {
    main <- full_gradient_run("main", retained_T, main_grid)
    warmup_seconds <- main$fit$stats$warmup_elapsed_time %||% NA_real_
  } else {
    bank_T <- pdmp_settings$anchor_bank_T
    warmup <- full_gradient_run("warmup", bank_T, grid_spacing = bank_T / 64)
    end <- materialize(warmup$fit)$terminal_states[[1L]]
    s <- subsampling(parameter_names)
    bank_seconds <- system.time(
      bank <- position_anchor_bank(stan_data, s$subsampling, parameter_names,
                                   warmup$online$grid_x,
                                   metric = s$bank_metric %||% "scaled"))[["elapsed"]]
    main <- do.call(run, c(list("main", main_grid), s$args, list(
      subsampling = s$subsampling, prior_standata = prior_standata,
      t0 = end$time, t_warmup = 0, T = end$time + retained_T,
      x0 = end$position, theta0 = end$physical_velocity,
      free0 = end$free_mask, stored_velocity0 = end$stored_frozen_velocity,
      flow_mean = end$flow$mean, flow_cov = end$flow$covariance,
      flow_specification = "frozen",
      # Subsampled bounds use finite-difference curvature (the package switches
      # to an exact Hessian-vector product where one exists).
      use_fd_hvp = TRUE, curvature_backend = "finite_difference",
      subsampling_anchor_selection_policy = "position",
      subsampling_position_anchor_bank_spec = bank,
      support_boundary = support_boundary_control(mode = "error"))))
    warmup_seconds <- (warmup$fit$stats$warmup_elapsed_time %||% NA_real_) +
      (warmup$fit$stats$main_elapsed_time %||% NA_real_) + bank_seconds
  }
  stats <- main$fit$stats
  # Numerical support-boundary recoveries (Bouncy Particle only; see the
  # supplement). They are not part of the PDMP, so every fit records how often
  # they happened, over the warmup stage as well for subsampled fits.
  recovery <- function(name) sum(vapply(
    list(main$fit, if (!is.null(subsampling)) warmup$fit),
    function(fit) if (is.null(fit)) 0 else fit$stats[[name]] %||% 0, numeric(1)))
  list(
    parameters = online_estimates(main$online, parameter_names, selectable),
    timing = c(total = main$elapsed + if (is.null(subsampling)) 0 else
                 warmup$elapsed + bank_seconds,
               warmup = warmup_seconds,
               sampling = stats$main_elapsed_time %||% NA_real_),
    counts = c(events = stats$main_events %||% NA_real_,
               gradient_calls = stats$main_gradient_calls %||% NA_real_,
               subset_evaluations =
                 stats$main_subsampling_subset_evaluations %||% NA_real_,
               boundary_recoveries = recovery("support_boundary_events"),
               boundary_refresh_failures =
                 recovery("support_boundary_refresh_failures")),
    # Everything below is stored so that any summary can be recomputed later
    # without rerunning: the raw online statistics of the retained phase (per
    # batch integrals and the state snapshots), every scalar sampler
    # statistic of each stage, the start, the adapted metric, and the anchors.
    online = main$online,
    stats = list(main = scalar_stats(main$fit$stats),
                 warmup = if (!is.null(subsampling))
                   scalar_stats(warmup$fit$stats)),
    start = list(x0 = x0, x0_sd = x0_sd),
    metric = terminal_metric(if (is.null(subsampling)) main$fit else warmup$fit),
    anchors = if (!is.null(subsampling)) bank$ordered_anchor_positions,
    # Timed batch integrals of the selectable coefficients for every phase in
    # chronological order (warmup, then for subsampled fits the anchor stage,
    # then the retained trajectory), for running estimates over computation
    # time (convergence_trace()).
    convergence = Filter(Negate(is.null), c(
      if (!is.null(subsampling)) list(
        warmup = timed_batches(warmup$online_warmup, parameter_names, selectable,
                               warmup$fit$stats$warmup_elapsed_time),
        anchor_stage = timed_batches(warmup$online, parameter_names, selectable,
          (warmup$fit$stats$main_elapsed_time %||% NA_real_) + bank_seconds)),
      list(warmup = if (is.null(subsampling))
             timed_batches(main$online_warmup, parameter_names, selectable,
                           main$fit$stats$warmup_elapsed_time),
           retained = timed_batches(main$online, parameter_names, selectable,
                                    stats$main_elapsed_time)))),
    anchor_bank_seconds = if (!is.null(subsampling)) bank_seconds,
    trajectory = if (!is.null(trajectory_dt))
      cbind(time = main$online$grid_times,
            `colnames<-`(t(main$online$grid_x[selectable, , drop = FALSE]),
                         parameter_names[selectable])))
}

# The batch integrals of the selectable coordinates of one phase and the
# duration of the phase in seconds (the sampler's phase timer).
timed_batches <- function(online, parameter_names, selectable, seconds) {
  if (is.null(online)) return(NULL)
  named <- function(m) `rownames<-`(m[selectable, , drop = FALSE],
                                    parameter_names[selectable])
  list(weight = online$batch_time, seconds = seconds %||% NA_real_,
       x = named(online$x), free = named(online$free))
}

# The same for a NIMBLE run (run_mcmc_batches()): per batch the number of
# iterations and the sums of the draws of the selectable coefficients
# (`columns`) and of their indicators, split into burn-in and retained
# batches; `seconds` is the run time of all iterations.
nimble_timed_batches <- function(batches, columns, nburnin, seconds,
                                 parameter_names = columns) {
  clean <- function(s) gsub("\\s+", "", s)
  cols <- match(clean(columns), clean(colnames(batches$sum)))
  stopifnot(!anyNA(cols))
  per_iteration <- seconds / sum(batches$iterations)
  phase <- ifelse(cumsum(batches$iterations) <= nburnin, "warmup", "retained")
  lapply(split(seq_along(phase), factor(phase, c("warmup", "retained"))),
         function(k) list(weight = batches$iterations[k],
                          seconds = per_iteration * sum(batches$iterations[k]),
                          x = `rownames<-`(t(batches$sum[k, cols, drop = FALSE]),
                                           parameter_names),
                          free = `rownames<-`(t(batches$nonzero[k, cols, drop = FALSE]),
                                              parameter_names)))
}

# Running estimates over computation time from the batches of a fit (its
# `convergence` element): after every batch, the posterior inclusion
# probabilities and posterior means of the selectable coefficients from all
# batches so far, and the computation time since the start of the run. Within
# a phase (PDMP warmup, anchor stage including the anchor bank, retained phase;
# NIMBLE burn-in and retained iterations) time is interpolated linearly in
# PDMP time or iterations from the duration of the phase, so it counts the
# phases themselves, as the run times do. With `from = "retained"` the
# estimates use the retained phase only, as the reported estimates do, while
# the time still counts from the start of the run.
# With `from = "second_half"` the estimate at time t uses only the second half
# of the computation time up to t (the first half is discarded as burn-in), so
# that an early transient does not stay in the estimates.
convergence_trace <- function(convergence,
                              from = c("start", "retained", "second_half")) {
  from <- match.arg(from)
  offset <- cumsum(c(0, vapply(convergence, `[[`, numeric(1), "seconds")))
  time <- unname(unlist(Map(function(p, o) o + p$seconds * cumsum(p$weight) / sum(p$weight),
                     convergence, offset[seq_along(convergence)])))
  phase <- rep(names(convergence), vapply(convergence, function(p) length(p$weight), 1L))
  phases <- if (from == "retained") convergence["retained"] else convergence
  weight <- cumsum(unlist(lapply(phases, `[[`, "weight")))
  # Start of the averaging window per batch: the first batch, or the first
  # batch that ends after half of the time so far.
  first <- if (from == "second_half")
    pmax(1L, findInterval(time / 2, time) + 1L) else rep(1L, length(weight))
  running <- function(field) {
    s <- do.call(cbind, lapply(phases, `[[`, field))
    cs <- cbind(0, t(apply(s, 1L, cumsum)))
    cw <- c(0, weight)
    k <- seq_along(weight)
    (cs[, k + 1L, drop = FALSE] - cs[, first, drop = FALSE]) /
      rep(cw[k + 1L] - cw[first], each = nrow(s))
  }
  list(time = if (from == "retained") time[phase == "retained"] else time,
       pip = running("free"), mean = running("x"))
}

# Time to convergence of Vogels et al. (2024, Equation 23): the first time at
# which a running accuracy measure is within `epsilon` of its final value,
#   T = min{t : |measure(t) - measure(final)| < epsilon},
# where `measure` is a vector over the times of convergence_trace().
time_to_convergence <- function(time, measure, epsilon) {
  hit <- which(abs(measure - measure[length(measure)]) < epsilon)
  if (length(hit)) time[hit[1L]] else NA_real_
}

# Area under the ROC curve of inclusion probabilities (ties count 1/2); NA
# when all coefficients are zero or all are nonzero.
auc <- function(pip, active) {
  if (all(active) || !any(active)) return(NA_real_)
  mean(outer(pip[active], pip[!active], ">") +
         0.5 * outer(pip[active], pip[!active], "=="))
}

# Times to convergence of one fit (Vogels et al., 2024): for the AUC of the
# running inclusion probabilities (within `epsilon_auc` of its final value) and
# for the root mean squared error of the running model-averaged posterior means
# (within the fraction `epsilon_rmse` of its final value). `truth` holds the
# true coefficients, named as the selectable coefficients.
convergence_times <- function(convergence, truth, epsilon_auc = 0.01,
                              epsilon_rmse = 0.05, from = "start") {
  tr <- convergence_trace(convergence, from)
  b <- truth[rownames(tr$pip)]
  stopifnot(!anyNA(b))
  auc_t <- apply(tr$pip, 2L, auc, active = b != 0)
  rmse_t <- sqrt(colMeans((tr$mean - b)^2))
  c(time_auc = if (all(is.na(auc_t))) NA_real_ else
      time_to_convergence(tr$time, auc_t, epsilon_auc),
    time_rmse = time_to_convergence(tr$time, rmse_t,
                                    epsilon_rmse * rmse_t[length(rmse_t)]),
    total_time = tr$time[length(tr$time)])
}

# Scalar entries of a sampler statistics list (drops long vectors).
scalar_stats <- function(stats) {
  keep <- vapply(stats, function(v) is.atomic(v) && length(v) == 1L, logical(1))
  unlist(stats[keep])
}

# Mean and diagonal covariance of the adapted flow at the end of a run.
terminal_metric <- function(fit) {
  flow <- tryCatch(materialize(fit)$terminal_states[[1L]]$flow,
                   error = function(e) NULL)
  if (is.null(flow)) return(NULL)
  list(mean = flow$mean, variance = diag(as.matrix(flow$covariance)))
}

# Posterior summaries from online statistics of the retained phase. By default
# the whole retained trajectory is used (the warmup precedes it, as the burn-in
# precedes NIMBLE's retained draws); `burnin` discards a leading fraction. The
# part used is split into `n_batches` equal batches for batch-means effective
# sample sizes.
#   mean    E[x], including the time at zero
#   pip     posterior inclusion probability, the fraction of time x != 0
# ESS is NA when a batch-means variance is zero or an edge never switched.
online_estimates <- function(online, parameter_names, selectable,
                             burnin = 0, n_batches = 40L) {
  keep <- (round(burnin * online$n_batches) + 1):online$n_batches
  group <- rep(seq_len(n_batches), each = length(keep) / n_batches)
  batch_sum <- function(m)
    t(apply(m[, keep, drop = FALSE], 1L, function(r) tapply(r, group, sum)))
  batch_time <- tapply(online$batch_time[keep], group, sum)
  total <- sum(batch_time)
  ess <- function(batch_means, variance) {
    v <- apply(batch_means, 1L, stats::var) / n_batches
    ifelse(v > 0, variance / v, NA_real_)
  }
  mean_batches <- sweep(batch_sum(online$x), 2L, batch_time, "/")
  pip_batches <- sweep(batch_sum(online$free), 2L, batch_time, "/")
  post_mean <- rowSums(batch_sum(online$x)) / total
  post_var <- pmax(0, rowSums(batch_sum(online$x2)) / total - post_mean^2)
  pip <- rowSums(batch_sum(online$free)) / total
  transitions <- rowSums(online$transitions[, keep, drop = FALSE])
  out <- data.frame(
    parameter = parameter_names, selectable = selectable, mean = post_mean,
    sd = sqrt(post_var), pip = ifelse(selectable, pip, NA_real_),
    ess_mean = ess(mean_batches, post_var),
    ess_pip = ifelse(selectable & transitions >= 2,
                     ess(pip_batches, pip * (1 - pip)), NA_real_),
    transitions = transitions)
  # Model size: the number of free selectable coordinates.
  n_fixed <- sum(!selectable)
  size_batches <- tapply(online$size[keep], group, sum) / batch_time - n_fixed
  size_mean <- sum(online$size[keep]) / total - n_fixed
  size_var <- max(0, sum(online$size2[keep]) / total - (size_mean + n_fixed)^2)
  rbind(out, data.frame(
    parameter = "model_size", selectable = FALSE, mean = size_mean,
    sd = sqrt(size_var), pip = NA_real_,
    ess_mean = ess(matrix(size_batches, 1L), size_var), ess_pip = NA_real_,
    transitions = NA_real_))
}

# Run fit_task(task) for every row of `tasks` on `workers` mirai daemons, in
# the order of `tasks`. Each result is saved as <out_dir>/fits/<task$id>.rds;
# tasks with an existing file are skipped, so an interrupted run can simply be
# restarted. Every fit runs in a fresh daemon that first sources `setup_file`
# (packages, Julia, fit_task) and exits after the fit, so memory does not build
# up across fits (compiled NIMBLE and Stan code, Julia's heap); the monitoring
# loop below starts a new daemon for each one that exits.
# Duration as "1h23m", "4m05s" or "12s".
format_seconds <- function(s) {
  if (!is.finite(s)) return("?")
  s <- round(s)
  if (s >= 3600) sprintf("%dh%02dm", s %/% 3600, s %% 3600 %/% 60) else
    if (s >= 60) sprintf("%dm%02ds", s %/% 60, s %% 60) else sprintf("%ds", s)
}

run_tasks <- function(tasks, setup_file, out_dir, workers) {
  fits_dir <- file.path(out_dir, "fits")
  dir.create(fits_dir, recursive = TRUE, showWarnings = FALSE)
  todo <- tasks[!file.exists(file.path(fits_dir, paste0(tasks$id, ".rds"))), ]
  message(nrow(todo), " of ", nrow(tasks), " fits to run on ", workers,
          " workers")
  if (!nrow(todo)) return(invisible(fits_dir))
  run_one <- function(task, fits_dir, dir, setup_file) {
    setwd(dir)
    Sys.setenv(SIMULATION_WORKER = "TRUE")
    source(setup_file)   # defines fit_task() in the daemon's global environment
    # mirai daemons use the L'Ecuyer-CMRG generator; use R's default, so that
    # every seed gives the same data and draws as in an ordinary R session.
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    result <- tryCatch(fit_task(task),
                       error = function(e) list(error = conditionMessage(e)))
    saveRDS(c(list(task = task), result),
            file.path(fits_dir, paste0(task$id, ".rds")))
    is.null(result$error)
  }
  environment(run_one) <- globalenv()   # resolve fit_task() on the daemon
  # One core per fit, so that wall times are comparable across methods and
  # parallel fits do not compete for cores: without these limits NIMBLE's
  # compiled code starts an OpenMP thread per core.
  Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1)
  mirai::daemons(workers, maxtasks = 1L)
  on.exit({
    # Every result is saved at this point. Julia's exit hooks can crash while R
    # shuts a worker down, so stop the workers directly.
    mirai::everywhere(tools::pskill(Sys.getpid(), tools::SIGKILL))
    mirai::daemons(0)
  })
  jobs <- mirai::mirai_map(split(todo, seq_len(nrow(todo))), run_one,
                           .args = list(fits_dir = fits_dir, dir = getwd(),
                                        setup_file = setup_file))
  # Progress over all fits of the study, including those finished by earlier
  # invocations, also written to <out_dir>/progress.txt for runs in the
  # background (e.g. `cat simulation_results/omrf_study/progress.txt`).
  finished_before <- nrow(tasks) - nrow(todo)
  started <- Sys.time()
  status <- ""
  shown <- -1L
  cli::cli_progress_bar(total = nrow(tasks), format = paste(
    "{cli::pb_bar} {cli::pb_current}/{cli::pb_total} fits | {status}"))
  repeat {
    done <- !vapply(jobs, mirai::unresolved, logical(1L))
    ok <- vapply(which(done), function(i) isTRUE(.subset2(jobs, i)$data),
                 logical(1L))
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    eta <- if (any(done)) elapsed / sum(done) * sum(!done) else NA
    status <- sprintf("%d failed | %s elapsed | ETA %s", sum(!ok),
                      format_seconds(elapsed), format_seconds(eta))
    # Redraw only when a fit finishes: in a log file every redraw is a line.
    if (sum(done) != shown) cli::cli_progress_update(
      set = finished_before + sum(done), force = TRUE)
    shown <- sum(done)
    writeLines(sprintf("%s  %d/%d fits done, %s", format(Sys.time()),
               finished_before + sum(done), nrow(tasks), status),
               file.path(out_dir, "progress.txt"))
    if (all(done)) break
    # Replace the daemons that exited after their fit, while fits are waiting.
    status_now <- mirai::status()
    n_new <- min(workers - status_now$connections,
                   status_now$mirai[["awaiting"]])
    if (n_new > 0) mirai::launch_local(n_new, maxtasks = 1L)
    Sys.sleep(10)
  }
  cli::cli_progress_done()
  message(sum(ok), " fits succeeded, ", sum(!ok), " failed")
  invisible(fits_dir)
}

# Collect saved fits into two tables: one row per fit, one row per parameter.
# Posterior summaries of a saved fit. PDMP summaries are recomputed from the
# stored batch integrals, so every fit uses the current online_estimates()
# whenever it was run; NIMBLE fits keep their stored summaries.
# Batch-means effective sample size of every column of a matrix of draws (the
# estimator used for the PDMP trajectories); NA for a constant column.
batch_means_ess <- function(draws, n_batches = 40L) {
  size <- nrow(draws) %/% n_batches
  draws <- draws[seq_len(size * n_batches), , drop = FALSE]
  group <- rep(seq_len(n_batches), each = size)
  means <- apply(draws, 2L, function(column) tapply(column, group, mean))
  batch_variance <- apply(means, 2L, stats::var) / n_batches
  ifelse(batch_variance > 0, apply(draws, 2L, stats::var) / batch_variance, NA_real_)
}

fit_parameters <- function(f) {
  p <- f$parameters
  if (is.null(f$online)) {
    # NIMBLE: effective sample size of the inclusion indicators from the draws.
    draws <- if (is.matrix(f$samples)) f$samples else f$samples$interactions
    if (!is.null(draws)) {
      ess <- batch_means_ess(draws != 0)
      m <- match(p$parameter, colnames(draws))
      p$ess_pip[!is.na(m)] <- ess[m[!is.na(m)]]
    }
    return(p)
  }
  keep <- p$parameter != "model_size"
  est <- online_estimates(f$online, p$parameter[keep], p$selectable[keep])
  extra <- setdiff(names(p), names(est))           # e.g. true_value
  est[extra] <- p[match(est$parameter, p$parameter), extra]
  est
}

# Phase timings of a saved fit. For subsampled fits the warmup is recomputed
# from the stored sampler timers (warmup stage: warmup plus anchor-stage
# sampling) plus the anchor-bank construction, so that fits run before
# fit_pdmp() timed them this way (when the warmup included Julia compilation)
# are summarised like the others; those lack the anchor-bank time (seconds).
fit_timing <- function(f) {
  timing <- f$timing %||% c(total = NA, warmup = NA, sampling = NA)
  w <- f$stats$warmup
  if (!is.null(w))
    timing[["warmup"]] <- as.numeric(w[["warmup_elapsed_time"]]) +
      as.numeric(w[["main_elapsed_time"]]) + (f$anchor_bank_seconds %||% 0)
  timing
}

collect_fits <- function(fits_dir) {
  # NIMBLE fits have no counts; fits from older script versions have fewer.
  count_template <- c(events = NA, gradient_calls = NA, subset_evaluations = NA,
                      boundary_recoveries = NA, boundary_refresh_failures = NA)
  # Fit files also hold raw online statistics and samples; read them one at a
  # time and keep only the summary tables, so memory use stays small.
  rows <- lapply(list.files(fits_dir, "\\.rds$", full.names = TRUE), function(file) {
    f <- readRDS(file)
    list(fit = data.frame(
           as.data.frame(f$task), error = f$error %||% NA_character_,
           t(fit_timing(f)),
           t(replace(count_template, names(f$counts), f$counts))),
         parameters = if (is.null(f$error)) cbind(id = f$task$id, fit_parameters(f)))
  })
  list(fits = do.call(rbind, lapply(rows, `[[`, "fit")),
       parameters = do.call(rbind, lapply(rows, `[[`, "parameters")))
}

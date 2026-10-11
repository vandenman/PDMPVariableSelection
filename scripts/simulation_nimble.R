# NIMBLE reversible-jump reference samplers for the simulation studies.
# Sourced on demand by simulation_omrf.R and simulation_clpm.R for NIMBLE tasks,
# so that PDMP-only workers never load NIMBLE.

require_nimble <- function() {
  if (!require("nimble", quietly = TRUE, character.only = TRUE)) {
    stop("Package 'nimble' is required.", call. = FALSE)
  }
  invisible(TRUE)
}

# Runs a compiled NIMBLE MCMC for `niter` iterations, burn-in included, and
# returns the retained draws (after `nburnin`), as runMCMC(nburnin = nburnin)
# would, together with the sums of all draws and of the indicators value != 0
# over about `n_batches` equal batches of the whole run (burn-in and retained
# iterations in separate batches). These give running estimates over the
# iterations, for the time to convergence.
nimble_batches <- 200L   # scripts/simulation_far_start.R uses finer batches
run_mcmc_batches <- function(cmcmc, niter, nburnin, seed, n_batches = nimble_batches) {
  draws <- as.matrix(nimble::runMCMC(cmcmc, niter = niter, nburnin = 0L,
    samplesAsCodaMCMC = FALSE, setSeed = seed))
  n_burnin <- max(1L, round(n_batches * nburnin / niter))
  batch <- c(cut(seq_len(nburnin), n_burnin, labels = FALSE),
             n_burnin + cut(seq_len(niter - nburnin), n_batches - n_burnin,
                            labels = FALSE))
  n_batches <- max(batch)
  list(samples = draws[-seq_len(nburnin), , drop = FALSE],
       batches = list(iterations = as.numeric(tabulate(batch, n_batches)),
                      sum = rowsum(draws, batch, reorder = TRUE),
                      nonzero = rowsum((draws != 0) + 0, batch, reorder = TRUE)))
}

clpm_nimble_code_generic <- function() {
  nimble::nimbleCode({
    for (j in 1:D) mu[j] ~ dnorm(0, sd = 1.5)
    # autoregressive effects N(0, 0.35^2); cross-lags N(0, slab sd^2), the slab
    # of the independent prior (the dependent-prior sampler replaces it)
    for (p in 1:P) for (q in 1:P) A[p, q] ~ dnorm(0, sd = A_prior_sd[p, q])
    for (p in 1:P) {
      sd_ri[p] ~ T(dnorm(0, sd = 1), 0, )
      sd_w1[p] ~ T(dnorm(0, sd = 1), 0, )
      sd_e[p] ~ T(dnorm(0, sd = 1), 0, )
    }
    Lcorr_ri[1:P, 1:P] ~ dlkj_corr_cholesky(eta = 2, p = P)
    Lcorr_w1[1:P, 1:P] ~ dlkj_corr_cholesky(eta = 2, p = P)
    Lcorr_e[1:P, 1:P] ~ dlkj_corr_cholesky(eta = 2, p = P)
    for (p in 1:P) for (q in 1:P) {
      U_ri[p, q] <- Lcorr_ri[p, q] * sd_ri[q]
      U_w1[p, q] <- Lcorr_w1[p, q] * sd_w1[q]
      U_e[p, q] <- Lcorr_e[p, q] * sd_e[q]
      Sigma_ri[p, q] <- inprod(U_ri[1:P, p], U_ri[1:P, q])
      Sigma_w1[p, q] <- inprod(U_w1[1:P, p], U_w1[1:P, q])
      Sigma_e[p, q] <- inprod(U_e[1:P, p], U_e[1:P, q])
      Cw[1, 1, p, q] <- Sigma_w1[p, q]
    }
    for (tt in 2:Tobs) {
      for (p in 1:P) for (q in 1:P)
        ACw[tt, p, q] <- inprod(A[p, 1:P], Cw[tt - 1, tt - 1, 1:P, q])
      for (p in 1:P) for (q in 1:P)
        Cw[tt, tt, p, q] <- inprod(ACw[tt, p, 1:P], A[q, 1:P]) + Sigma_e[p, q]
      for (ss in 1:(tt - 1)) for (p in 1:P) for (q in 1:P) {
        Cw[tt, ss, p, q] <- inprod(A[p, 1:P], Cw[tt - 1, ss, 1:P, q])
        Cw[ss, tt, p, q] <- Cw[tt, ss, q, p]
      }
    }
    for (tt in 1:Tobs) for (ss in 1:Tobs) for (p in 1:P) for (q in 1:P)
      Sigma_z[(tt - 1) * P + p, (ss - 1) * P + q] <-
        Sigma_ri[p, q] + Cw[tt, ss, p, q] +
        equals((tt - 1) * P + p, (ss - 1) * P + q) * jitter
    Omega_z[1:D, 1:D] <- inverse(Sigma_z[1:D, 1:D])
    for (n in 1:N) z[n, 1:D] ~ dmnorm(mean = mu[1:D], prec = Omega_z[1:D, 1:D])
  })
}

prepare_clpm_data_generic <- function(z, P, Tobs, jitter = 1e-8) {
  z <- as.matrix(z)
  N <- nrow(z); D <- ncol(z)
  stopifnot(N >= 1L, P >= 1L, Tobs >= 2L, D == P * Tobs, jitter >= 0)
  list(N = as.integer(N), P = as.integer(P), Tobs = as.integer(Tobs),
       D = as.integer(D), z = z, jitter = as.numeric(jitter))
}

clpm_nimble_inits_generic <- function(P, Tobs) {
  list(mu = numeric(P * Tobs), A = diag(0.30, P), sd_ri = rep(0.60, P),
       sd_w1 = rep(0.70, P), sd_e = rep(0.50, P), Lcorr_ri = diag(P),
       Lcorr_w1 = diag(P), Lcorr_e = diag(P))
}

# Two implementations of the RI-CLPM likelihood with identical posteriors:
#   "person"      the natural implementation given the model definition: every
#                 person's pT responses are multivariate normal, one dmnorm node
#                 per person (cost O(N) per evaluation). Kept to validate the
#                 optimized version (compare_clpm_nimble_likelihoods()).
#   "sufficient"  the optimized version used in the study: one node whose
#                 density uses the sample mean zbar and the centred scatter
#                 matrix S, computed once (cost independent of N):
#                   log L = N/2 log|Omega| - N D/2 log(2 pi)
#                           - 1/2 [tr(Omega S) + N (zbar - mu)' Omega (zbar - mu)],
#                 with the same precision matrix Omega_z (jitter included).
dmnorm_suff <- if (requireNamespace("nimble", quietly = TRUE)) {
  nimble::nimbleFunction(
    run = function(x = double(0), mean = double(1), prec = double(2),
                   N = double(0), zbar = double(1), S = double(2),
                   log = integer(0, default = 0)) {
      returnType(double(0))
      D <- length(mean)
      U <- chol(prec)
      log_det_prec <- 2 * sum(log(diag(U)))
      diff <- zbar - mean
      quad <- sum(prec * S) + N * sum(diff * (prec %*% diff)[, 1])
      ll <- 0.5 * N * log_det_prec - 0.5 * N * D * log(2 * pi) - 0.5 * quad
      if (log) return(ll)
      return(exp(ll))
    })
} else NULL
rmnorm_suff <- if (requireNamespace("nimble", quietly = TRUE)) {
  nimble::nimbleFunction(
    run = function(n = integer(0), mean = double(1), prec = double(2),
                   N = double(0), zbar = double(1), S = double(2)) {
      returnType(double(0))
      return(0)    # the node is data; it is never simulated
    })
} else NULL

register_clpm_nimble_distributions <- function() {
  nimble::registerDistributions(list(dmnorm_suff = list(
    BUGSdist = "dmnorm_suff(mean, prec, N, zbar, S)",
    types = c("value = double(0)", "mean = double(1)", "prec = double(2)",
              "N = double(0)", "zbar = double(1)", "S = double(2)"))),
    verbose = FALSE)
  invisible(TRUE)
}

# The model with its last statement (one dmnorm node per person) replaced by
# the sufficient-statistic node.
clpm_nimble_code_sufficient <- function() {
  code <- clpm_nimble_code_generic()
  last <- length(code)
  stopifnot(grepl("dmnorm", deparse(code[[last]])[1]))
  code[[last]] <- quote(zsuff ~ dmnorm_suff(mean = mu[1:D], prec = Omega_z[1:D, 1:D],
                                            N = N, zbar = zbar[1:D], S = S[1:D, 1:D]))
  code
}

build_clpm_nimble_model_generic <- function(clpm_data, inits = NULL,
                                            likelihood = c("sufficient", "person"),
                                            cross_lag_sd = 0.35) {
  likelihood <- match.arg(likelihood)
  inits <- inits %||% clpm_nimble_inits_generic(clpm_data$P, clpm_data$Tobs)
  constants <- clpm_data[c("N", "P", "Tobs", "D", "jitter")]
  constants$A_prior_sd <- matrix(cross_lag_sd, clpm_data$P, clpm_data$P)
  diag(constants$A_prior_sd) <- 0.35
  if (likelihood == "person")
    return(nimble::nimbleModel(code = clpm_nimble_code_generic(),
      constants = constants, data = clpm_data["z"], inits = inits,
      calculate = FALSE))
  register_clpm_nimble_distributions()
  zbar <- colMeans(clpm_data$z)
  nimble::nimbleModel(code = clpm_nimble_code_sufficient(),
    constants = c(constants, list(zbar = zbar,
                                  S = crossprod(sweep(clpm_data$z, 2, zbar)))),
    data = list(zsuff = 0), inits = inits, calculate = FALSE)
}

# Validation: the two likelihoods give the same log posterior (up to rounding)
# at random parameter values. Returns the largest absolute difference.
compare_clpm_nimble_likelihoods <- function(clpm_data, n_values = 20L, seed = 1L) {
  person <- build_clpm_nimble_model_generic(clpm_data, likelihood = "person")
  sufficient <- build_clpm_nimble_model_generic(clpm_data, likelihood = "sufficient")
  set.seed(seed)
  diffs <- vapply(seq_len(n_values), function(r) {
    values <- list(mu = stats::rnorm(clpm_data$D, sd = 0.2),
                   A = diag(0.35, clpm_data$P) +
                     matrix(stats::rnorm(clpm_data$P^2, sd = 0.1), clpm_data$P),
                   sd_ri = stats::runif(clpm_data$P, 0.4, 0.8),
                   sd_w1 = stats::runif(clpm_data$P, 0.5, 0.9),
                   sd_e = stats::runif(clpm_data$P, 0.4, 0.6))
    for (name in names(values)) {
      person[[name]] <- values[[name]]
      sufficient[[name]] <- values[[name]]
    }
    sufficient$calculate() - person$calculate()
  }, numeric(1))
  max(abs(diffs))
}

# The dependent CLPM sampler uses the same fixed-dimensional reversible-jump
# correction as the previous implementation. Each call sweeps over all target
# coefficients (like configureRJ does for the independent prior, which gives
# every coefficient its own sampler), updating them one at a time.  It is defined only when the
# NIMBLE module is loaded, never when a PDMP child sources simulation_pdmp.R.
sampler_clpm_dependent_rj <- if (requireNamespace("nimble", quietly = TRUE)) {
  nimble::nimbleFunction(
    contains = nimble::sampler_BASE,
    setup = function(model, mvSaved, target, control) {
      targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
      cc <- nimble:::mcmc_determineCalcAndCopyNodes(model, targetAsScalar)
      calcNodesNoSelf <- cc$calcNodesNoSelf
      copyNodesDeterm <- cc$copyNodesDeterm
      copyNodesStoch <- cc$copyNodesStoch
      proposalScale <- extractControlElement(control, "proposalScale", 0.35)
      rwScale <- extractControlElement(control, "rwScale", proposalScale)
      slabSD <- extractControlElement(control, "slabSD", 0.35)
      slabCorrelation <- extractControlElement(control, "slabCorrelation", 0.25)
      betaA <- extractControlElement(control, "betaA", 1)
      betaB <- extractControlElement(control, "betaB", 4)
      baseSD <- extractControlElement(control, "baseSD", 0.35)
      zeroTol <- extractControlElement(control, "zeroTol", 1e-12)
      activeUpdateProbability <- extractControlElement(control, "activeUpdateProbability", 0.5)
      d <- length(targetAsScalar)
      u <- slabSD * slabSD * (1 - slabCorrelation)
      v <- slabSD * slabSD * slabCorrelation
      if (d < 1 || proposalScale <= 0 || rwScale <= 0 || slabSD <= 0 ||
          slabCorrelation < 0 || slabCorrelation >= 1 || activeUpdateProbability < 0 ||
          activeUpdateProbability >= 1 || u <= 0 || u + d * v <= 0)
        stop("Invalid dependent CLPM NIMBLE sampler controls")
      timesRan <- 0; timesAccepted <- 0
    },
    run = function() {
      proposalValue <- numeric(1)
      for (i in 1:d) {
        oldValue <- values(model, targetAsScalar[i])[1]
        oldCorrection <- logPriorCorrection(); oldActive <- abs(oldValue) > zeroTol
        logProposalRatio <- 0
        if (oldActive) {
          if (runif(1) < activeUpdateProbability) newValue <- rnorm(1, oldValue, rwScale)
          else { newValue <- 0; logProposalRatio <- dnorm(oldValue, 0, proposalScale, log = TRUE) - log(1 - activeUpdateProbability) }
        } else {
          newValue <- rnorm(1, 0, proposalScale)
          logProposalRatio <- log(1 - activeUpdateProbability) - dnorm(newValue, 0, proposalScale, log = TRUE)
        }
        proposalValue[1] <- newValue; values(model, targetAsScalar[i]) <<- proposalValue
        logMHR <- model$calculateDiff(targetAsScalar[i])
        jump <- FALSE
        # NIMBLE cannot compile the short-circuit `&&` expression that was
        # previously used for this acceptance rule.  Nested scalar comparisons
        # preserve the original finite-only guard: NaN and both infinities
        # reject, and only then is the log-uniform comparison evaluated.
        logMHR <- logMHR + model$calculateDiff(calcNodesNoSelf) +
          logPriorCorrection() - oldCorrection + logProposalRatio
        if (logMHR > -Inf) {
          if (logMHR < Inf) {
            logUniform <- log(runif(1))
            if (logUniform < logMHR) jump <- TRUE
          }
        }
        if (jump) {
          nimCopy(from = model, to = mvSaved, row = 1, nodes = targetAsScalar, logProb = TRUE)
          nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
          nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
        } else {
          nimCopy(from = mvSaved, to = model, row = 1, nodes = targetAsScalar, logProb = TRUE)
          nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
          nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
        }
        timesRan <<- timesRan + 1; if (jump) timesAccepted <<- timesAccepted + 1
      }
    },
    methods = list(
      logPriorCorrection = function() {
        returnType(double())
        k <- 0; sumBeta <- 0; sumBetaSq <- 0; baseLogDensity <- 0
        for (j in 1:d) {
          beta <- values(model, targetAsScalar[j])[1]
          baseLogDensity <- baseLogDensity + dnorm(beta, 0, baseSD, log = TRUE)
          if (abs(beta) > zeroTol) { k <- k + 1; sumBeta <- sumBeta + beta; sumBetaSq <- sumBetaSq + beta * beta }
        }
        modelLogPrior <- lgamma(betaA + k) + lgamma(betaB + d - k) - lgamma(betaA + betaB + d)
        slabLogDensity <- 0
        if (k > 0) {
          # Keep this as ordinary scalar control flow.  NIMBLE 1.4.x can
          # misclassify the literal branch of an inline `if` expression when
          # it is nested in a custom sampler method.
          logDet <- 0
          if (k > 1) logDet <- (k - 1) * log(u)
          logDet <- logDet + log(u + k * v)
          quad <- sumBetaSq / u - (v * sumBeta * sumBeta) / (u * (u + k * v))
          slabLogDensity <- -0.5 * (k * log(2 * pi) + logDet + quad)
        }
        return(modelLogPrior + slabLogDensity - baseLogDensity)
      },
      # NIMBLE 1.4.2 cannot lower `NaN` as a scalar return expression in this
      # custom sampler.  Zero is the natural empty-sample convention and the
      # method is only used for reporting after at least one update in normal
      # runs.
      getAcceptanceRate = function() { returnType(double()); if (timesRan == 0) return(0.0); return(timesAccepted / timesRan) },
      reset = function() { timesRan <<- 0; timesAccepted <<- 0 }
    )
  )
} else NULL

configure_clpm_nimble_generic <- function(model, target_nodes = character(), prior_prob = 0.20,
    proposal_scale = 0.35, use_conjugacy = FALSE, prior_type = c("independent", "dependent"),
    beta_binomial_a = 1, beta_binomial_b = 4, slab_sd = 0.35, slab_correlation = 0.25,
    active_rw_scale = proposal_scale, active_update_probability = 0.5) {
  prior_type <- match.arg(prior_type)
  conf <- nimble::configureMCMC(model, monitors = c("A", "mu", "sd_ri", "sd_w1", "sd_e", "Lcorr_ri", "Lcorr_w1", "Lcorr_e"),
    enableWAIC = FALSE, print = FALSE, useConjugacy = use_conjugacy)
  if (length(target_nodes)) {
    if (identical(prior_type, "dependent")) {
      conf$removeSamplers(target_nodes)
      conf$addSampler(target = target_nodes, type = sampler_clpm_dependent_rj,
        control = list(proposalScale = proposal_scale, rwScale = active_rw_scale,
          slabSD = slab_sd, slabCorrelation = slab_correlation, betaA = beta_binomial_a,
          betaB = beta_binomial_b, baseSD = slab_sd, activeUpdateProbability = active_update_probability),
        name = "clpm_dependent_rj", print = FALSE)
    } else nimble::configureRJ(conf, targetNodes = target_nodes, priorProb = prior_prob,
      control = list(mean = 0, scale = proposal_scale))
  }
  conf
}

clpm_off_diagonal_nodes <- function(P) unlist(lapply(seq_len(P), function(p)
  sprintf("A[%d,%d]", p, setdiff(seq_len(P), p))), use.names = FALSE)

normalize_nimble_node_names <- function(x) gsub("\\s+", "", x)

run_clpm_nimble_generic <- function(
    clpm_data,
    niter = 4000L,
    nburnin = 1000L,
    thin = 1L,
    seed = 20260702L,
    target_nodes = clpm_off_diagonal_nodes(clpm_data$P),
    prior_prob = 0.20,
    proposal_scale = 0.35,
    prior_type = c("independent", "dependent"),
    beta_binomial_a = 1,
    beta_binomial_b = 4,
    slab_sd = 0.35,
    slab_correlation = 0.25,
    active_rw_scale = proposal_scale,
    active_update_probability = 0.5,
    compile = TRUE,
    likelihood = c("sufficient", "person"),
    inits = NULL) {
  require_nimble()
  prior_type <- match.arg(prior_type)
  likelihood <- match.arg(likelihood)
  stopifnot(requireNamespace("coda", quietly = TRUE))
  model <- build_clpm_nimble_model_generic(clpm_data, inits = inits,
                                           likelihood = likelihood,
                                           cross_lag_sd = slab_sd)
  initial_log_prob <- model$calculate()
  conf <- configure_clpm_nimble_generic(
    model,
    target_nodes = target_nodes,
    prior_prob = prior_prob,
    proposal_scale = proposal_scale,
    prior_type = prior_type,
    beta_binomial_a = beta_binomial_a,
    beta_binomial_b = beta_binomial_b,
    slab_sd = slab_sd,
    slab_correlation = slab_correlation,
    active_rw_scale = active_rw_scale,
    active_update_probability = active_update_probability
  )
  mcmc <- nimble::buildMCMC(conf)

  compile_sec <- 0
  if (isTRUE(compile)) {
    compiled <- NULL
    compile_sec <- system.time({
      compiled <- nimble::compileNimble(model, mcmc)
    })[["elapsed"]]
    run_mcmc <- compiled$mcmc
  } else {
    run_mcmc <- mcmc
  }

  stopifnot(as.integer(thin) == 1L)
  sample_sec <- system.time({
    run <- run_mcmc_batches(run_mcmc, as.integer(niter), as.integer(nburnin),
                            as.integer(seed))
  })[["elapsed"]]
  samples <- run$samples
  if (isTRUE(compile)) {
    try(nimble::clearCompiled(compiled), silent = TRUE)
    rm(run_mcmc, compiled, mcmc, conf, model)
    invisible(gc())
  }
  normalized_columns <- normalize_nimble_node_names(colnames(samples))
  normalized_targets <- normalize_nimble_node_names(target_nodes)
  target_indices <- match(normalized_targets, normalized_columns)
  if (anyNA(target_indices)) {
    stop("Missing monitored target nodes: ",
         paste(target_nodes[is.na(target_indices)], collapse = ", "), call. = FALSE)
  }
  target_samples <- samples[, target_indices, drop = FALSE]
  colnames(target_samples) <- normalized_targets

  list(
    model = if (isTRUE(compile)) NULL else model,
    samples = samples,
    target_samples = target_samples,
    batches = run$batches,
    initial_log_prob = initial_log_prob,
    compile_sec = compile_sec,
    sample_sec = sample_sec,
    target_ess = as.numeric(coda::effectiveSize(target_samples)),
    target_inclusion = colMeans(target_samples != 0),
    target_transitions = apply(target_samples != 0, 2L, function(x) sum(diff(x) != 0)),
    control = list(
      likelihood = likelihood,
      niter = as.integer(niter),
      nburnin = as.integer(nburnin),
      thin = as.integer(thin),
      seed = as.integer(seed),
      target_nodes = normalized_targets,
      prior_type = prior_type,
      prior_prob = prior_prob,
      proposal_scale = proposal_scale,
      beta_binomial_a = beta_binomial_a,
      beta_binomial_b = beta_binomial_b,
      slab_sd = slab_sd,
      slab_correlation = slab_correlation,
      active_rw_scale = active_rw_scale,
      active_update_probability = active_update_probability
    )
  )
}
run_clpm_nimble_rj <- function(stan_data, truth = NULL, niter = 4000L, nburnin = 1000L,
                               thin = 1L, nchains = 1L, seed = 20260702L,
	                               target_nodes = clpm_off_diagonal_nodes(stan_data$P),
	                               prior_prob = 0.20, proposal_scale = 0.35,
	                               prior_type = c("independent", "dependent"),
	                               beta_binomial_a = 1, beta_binomial_b = 4,
	                               slab_sd = 0.35, slab_correlation = 0.25,
	                               active_rw_scale = proposal_scale,
	                               active_update_probability = 0.5,
	                               compile = TRUE,
	                               likelihood = c("sufficient", "person"),
	                               inits = NULL) {
  prior_type <- match.arg(prior_type)
  likelihood <- match.arg(likelihood)
  if (as.integer(nchains) != 1L) {
    stop("The generic CLPM NIMBLE compatibility wrapper currently supports one chain.", call. = FALSE)
  }
  clpm_data <- prepare_clpm_data_generic(
    z = stan_data$z,
    P = as.integer(stan_data$P),
    Tobs = as.integer(stan_data$T)
  )
  out <- run_clpm_nimble_generic(
    clpm_data = clpm_data,
    niter = niter,
    nburnin = nburnin,
    thin = thin,
    seed = seed,
    target_nodes = target_nodes,
    prior_prob = prior_prob,
    proposal_scale = proposal_scale,
    prior_type = prior_type,
    beta_binomial_a = beta_binomial_a,
    beta_binomial_b = beta_binomial_b,
    slab_sd = slab_sd,
    slab_correlation = slab_correlation,
    active_rw_scale = active_rw_scale,
    active_update_probability = active_update_probability,
    compile = compile,
    likelihood = likelihood,
    inits = inits
  )

  list(
    status = "fit",
    method = "nimble_rj",
    prior_type = prior_type,
    truth = truth,
    control = c(
      out$control,
      list(
        post_warmup = max(0L, as.integer(niter) - as.integer(nburnin)),
        nchains = as.integer(nchains),
        compile = compile
      )
    ),
    fit = structure(
      list(
        stats = list(
          elapsed_time = out$sample_sec,
          compile_time = out$compile_sec,
          ct_ess = matrix(as.numeric(out$target_ess), nrow = 1L)
        ),
        samples = out$samples,
        target_samples = out$target_samples,
        target_transitions = out$target_transitions,
        batches = out$batches,
        posterior_summary = data.frame(
          parameter = colnames(out$target_samples),
          posterior_mean = colMeans(out$target_samples),
          cri_low = apply(out$target_samples, 2, stats::quantile, probs = 0.025),
          cri_high = apply(out$target_samples, 2, stats::quantile, probs = 0.975),
          inclusion_probability = out$target_inclusion,
          ess = as.numeric(out$target_ess),
          stringsAsFactors = FALSE
        ),
        model = out$model
      ),
      class = "nimble_clpm_rj_result"
    ),
    unc_names = colnames(out$target_samples),
    can_stick = rep(TRUE, ncol(out$target_samples))
  )
}

summarise_clpm_nimble_rj <- function(result) {
  summary <- result$fit$posterior_summary
  summary$elapsed_time <- result$fit$stats$elapsed_time
  summary$ess_per_sec <- summary$ess / result$fit$stats$elapsed_time
  summary
}


# OMRF NIMBLE reference sampler helpers ----
require_omrf_nimble <- function() {
  if (!require("nimble", quietly = TRUE, character.only = TRUE)) {
    stop("Package 'nimble' is required for the OMRF RJMCMC reference.", call. = FALSE)
  }
  invisible(TRUE)
}


dlogitbeta <- if (requireNamespace("nimble", quietly = TRUE)) {
  nimble::nimbleFunction(
  run = function(x = double(0), alpha = double(0), beta = double(0),
                 log = integer(0, default = 0)) {
    returnType(double(0))
    # The omitted beta-function term is constant in x and cancels in MCMC ratios.
    log_density <- alpha * x - (alpha + beta) * log1p(exp(x))
    if (log) return(log_density)
    return(exp(log_density))
  }
  )
} else {
  NULL
}

rlogitbeta <- if (requireNamespace("nimble", quietly = TRUE)) {
  nimble::nimbleFunction(
  run = function(n = integer(0), alpha = double(0), beta = double(0)) {
    returnType(double(0))
    probability <- rbeta(1, alpha, beta)
    return(log(probability / (1.0 - probability)))
  }
  )
} else {
  NULL
}

register_omrf_nimble_distributions <- function() {
  nimble::registerDistributions(list(
    dlogitbeta = list(
      BUGSdist = "dlogitbeta(alpha, beta)",
      Rdist = "dlogitbeta(alpha, beta)",
      types = c("value = double(0)", "alpha = double(0)", "beta = double(0)")
    )
  ))
  invisible(TRUE)
}

# Two implementations of the OMRF pseudolikelihood with identical targets:
#   "dense"     the natural implementation given the model definition: each
#               person-node field is the inner product of all E interactions
#               with that node's edge covariates, and the pseudolikelihood
#               enters per person. NIMBLE then recomputes all N P fields (at
#               O(E) each) and all N person totals whenever any single edge is
#               updated. Kept to validate the optimized version
#               (compare_omrf_nimble_implementations()).
#   "incident"  the optimized version used in the study: each field sums only
#               the P - 1 edges incident to its node, and the pseudolikelihood
#               enters per person-node, so an edge update recomputes 2N fields
#               and person-node terms. The chain is the same, about 4 times
#               faster at P = 10.
omrf_nimble_likelihood_code <- function(implementation) {
  switch(match.arg(implementation, c("incident", "dense")),
    dense = quote(for (n in 1:N) {
      for (node in 1:P) {
        field[n, node] <- inprod(interactions[1:E], edge_covariate[n, node, 1:E])
        for (u in 1:(K - 1)) {
          denominator_term[n, node, u] <- exp(thresholds[node, u] + u * field[n, node])
        }
        denominator[n, node] <- 1 + sum(denominator_term[n, node, 1:(K - 1)])
        log_pseudolik_node[n, node] <-
          threshold_expanded[node, X[n, node] + 1] +
          X[n, node] * field[n, node] - log(denominator[n, node])
      }
      negative_log_pseudolik[n] <- -sum(log_pseudolik_node[n, 1:P])
      zeros[n] ~ dpois(negative_log_pseudolik[n])
    }),
    incident = quote(for (n in 1:N) {
      for (node in 1:P) {
        # NIMBLE does not index by a vector of indices: one term per incident edge
        for (k in 1:Q) {
          incident_term[n, node, k] <-
            interactions[incident_edge[node, k]] * neighbour_x[n, node, k]
        }
        field[n, node] <- sum(incident_term[n, node, 1:Q])
        for (u in 1:(K - 1)) {
          denominator_term[n, node, u] <- exp(thresholds[node, u] + u * field[n, node])
        }
        denominator[n, node] <- 1 + sum(denominator_term[n, node, 1:(K - 1)])
        log_pseudolik_node[n, node] <-
          threshold_expanded[node, X[n, node] + 1] +
          X[n, node] * field[n, node] - log(denominator[n, node])
        # Poisson-zeros trick: the log density of zeros[n, node] is log_pseudolik_node
        zeros[n, node] ~ dpois(-log_pseudolik_node[n, node])
      }
    }))
}

omrf_nimble_code <- function(implementation = c("incident", "dense")) {
  register_omrf_nimble_distributions()
  code <- nimble::nimbleCode({
    rho ~ dbeta(beta_binomial_a, beta_binomial_b)

    for (node in 1:P) {
      threshold_expanded[node, 1] <- 0
      for (u in 1:(K - 1)) {
        thresholds[node, u] ~ dlogitbeta(prior_threshold_alpha, prior_threshold_beta)
        threshold_expanded[node, u + 1] <- thresholds[node, u]
      }
    }

    for (edge in 1:E) {
      gamma[edge] ~ dbern(rho)
      slab_z[edge] ~ dt(0, sigma = 1, df = 1)
      interactions[edge] <- gamma[edge] * prior_cauchy_scale * slab_z[edge]
    }

    LIKELIHOOD
  })
  code[[length(code)]] <- omrf_nimble_likelihood_code(implementation)
  code
}

omrf_nimble_code_node <- function(implementation = c("incident", "dense")) {
  register_omrf_nimble_distributions()
  code <- nimble::nimbleCode({
    rho ~ dbeta(beta_binomial_a, beta_binomial_b)

    for (node in 1:P) {
      node_scale[node] ~ dlnorm(0, sdlog = prior_node_scale_sd)
      threshold_expanded[node, 1] <- 0
      for (u in 1:(K - 1)) {
        thresholds[node, u] ~ dlogitbeta(prior_threshold_alpha, prior_threshold_beta)
        threshold_expanded[node, u + 1] <- thresholds[node, u]
      }
    }

    for (edge in 1:E) {
      gamma[edge] ~ dbern(rho)
      edge_sd[edge] <- slab_base_scale *
        sqrt(node_scale[edge_i[edge]] * node_scale[edge_j[edge]])
      slab_z[edge] ~ dnorm(0, sd = 1)
      interactions[edge] <- gamma[edge] * edge_sd[edge] * slab_z[edge]
    }

    LIKELIHOOD
  })
  code[[length(code)]] <- omrf_nimble_likelihood_code(implementation)
  code
}

omrf_edge_covariate <- function(X, edges) {
  N <- nrow(X)
  P <- ncol(X)
  E <- nrow(edges)
  out <- array(0, dim = c(N, P, E))
  for (edge in seq_len(E)) {
    i <- edges$i[[edge]]
    j <- edges$j[[edge]]
    out[, i, edge] <- X[, j]
    out[, j, edge] <- X[, i]
  }
  out
}

omrf_nimble_inits <- function(stan_data, prior_type = "beta_binomial_independent",
                              seed = 20260702L, start_label = NULL) {
  set.seed(seed)
  E <- stan_data$P * (stan_data$P - 1L) / 2L
  gamma <- if (is.null(start_label)) {
    stats::rbinom(E, 1L, 0.2)
  } else {
    start_label <- match.arg(start_label,
      c("empty", "sparse", "random", "complete"))
    switch(start_label,
      empty = integer(E),
      sparse = replace(integer(E), seq_len(max(1L, ceiling(.10 * E))), 1L),
      random = replace(integer(E), sample.int(E, ceiling(.50 * E)), 1L),
      complete = rep(1L, E))
  }
  inits <- list(
    rho = 0.2,
    gamma = gamma,
    thresholds = matrix(0, stan_data$P, stan_data$K - 1L)
  )
  if (identical(prior_type, "beta_binomial_node")) {
    inits$node_scale <- rep(1, stan_data$P)
    # Fixed-dimensional augmentation for the stochastic node-specific slab.
    # Marginally, gamma * edge_sd * slab_z is the intended point-mass/normal
    # mixture, while slab_z retains its standard-normal pseudo-prior when the
    # edge is inactive.
    inits$slab_z <- stats::rnorm(E)
  } else {
    # A moderate valid initial value is preferable to drawing directly from
    # the heavy-tailed Cauchy pseudo-prior.
    inits$slab_z <- stats::rnorm(E, 0, 0.10)
  }
  inits
}

# Incident edges of every node and the corresponding neighbour responses, for
# the "incident" implementation.
omrf_incident_edges <- function(X, edges) {
  P <- ncol(X)
  incident <- t(vapply(seq_len(P), function(node)
    which(edges$i == node | edges$j == node), integer(P - 1)))
  neighbour <- t(vapply(seq_len(P), function(node) {
    e <- incident[node, ]
    ifelse(edges$i[e] == node, edges$j[e], edges$i[e])
  }, integer(P - 1)))
  neighbour_x <- array(0, c(nrow(X), P, P - 1))
  for (node in seq_len(P)) neighbour_x[, node, ] <- X[, neighbour[node, ]]
  list(incident_edge = incident, neighbour_x = neighbour_x, Q = as.integer(P - 1))
}

omrf_nimble_model <- function(stan_data, beta_binomial_a = 1, beta_binomial_b = 4,
                              prior_type = "beta_binomial_independent",
                              inits = NULL,
                              implementation = c("incident", "dense")) {
  implementation <- match.arg(implementation)
  edges <- edge_index(stan_data$P)
  constants <- list(
    N = as.integer(stan_data$N),
    P = as.integer(stan_data$P),
    K = as.integer(stan_data$K),
    E = as.integer(nrow(edges)),
    edge_i = as.integer(edges$i),
    edge_j = as.integer(edges$j),
    X = stan_data$X,
    beta_binomial_a = beta_binomial_a,
    beta_binomial_b = beta_binomial_b,
    prior_cauchy_scale = stan_data$prior_cauchy_scale,
    prior_node_scale_sd = stan_data$prior_node_scale_sd,
    slab_base_scale = stan_data$slab_base_scale,
    prior_threshold_alpha = stan_data$prior_threshold_alpha,
    prior_threshold_beta = stan_data$prior_threshold_beta
  )
  if (implementation == "dense") {
    constants$edge_covariate <- omrf_edge_covariate(stan_data$X, edges)
    data <- list(zeros = rep(0L, stan_data$N))
  } else {
    constants <- c(constants, omrf_incident_edges(stan_data$X, edges))
    data <- list(zeros = matrix(0L, stan_data$N, stan_data$P))
  }
  inits <- inits %||% omrf_nimble_inits(stan_data, prior_type)
  code <- if (identical(prior_type, "beta_binomial_node")) {
    omrf_nimble_code_node(implementation)
  } else if (identical(prior_type, "beta_binomial_independent")) {
    omrf_nimble_code(implementation)
  } else {
    stop("Unknown OMRF NIMBLE prior type: ", prior_type, call. = FALSE)
  }
  nimble::nimbleModel(
    code = code,
    constants = constants,
    data = data,
    inits = inits,
    calculate = FALSE
  )
}

# Validation: the dense and incident implementations give the same log posterior
# (up to rounding) at random parameter values. Returns the largest absolute
# difference.
compare_omrf_nimble_implementations <- function(stan_data,
    prior_type = "beta_binomial_independent", n_values = 10L, seed = 1L) {
  inits <- omrf_nimble_inits(stan_data, prior_type, seed)
  dense <- omrf_nimble_model(stan_data, prior_type = prior_type, inits = inits,
                             implementation = "dense")
  incident <- omrf_nimble_model(stan_data, prior_type = prior_type, inits = inits,
                                implementation = "incident")
  E <- stan_data$P * (stan_data$P - 1) / 2
  set.seed(seed)
  diffs <- vapply(seq_len(n_values), function(r) {
    values <- list(rho = stats::runif(1, 0.1, 0.5), gamma = stats::rbinom(E, 1, 0.4),
                   slab_z = stats::rnorm(E, 0, 0.3),
                   thresholds = matrix(stats::rnorm(stan_data$P * (stan_data$K - 1)), stan_data$P))
    if (identical(prior_type, "beta_binomial_node"))
      values$node_scale <- stats::rlnorm(stan_data$P, 0, 0.3)
    for (name in names(values)) {
      dense[[name]] <- values[[name]]
      incident[[name]] <- values[[name]]
    }
    dense$calculate() - incident$calculate()
  }, numeric(1))
  max(abs(diffs))
}

# Latent slab value of an edge: while the edge is excluded (gamma = 0) the
# value does not enter the likelihood, so its full conditional is its prior and
# it is drawn exactly from it (no likelihood evaluation); while the edge is
# included it gets the adaptive random-walk update.
sampler_omrf_slab_z <- if (requireNamespace("nimble", quietly = TRUE)) {
  nimble::nimbleFunction(
    name = "sampler_omrf_slab_z", contains = nimble::sampler_BASE,
    setup = function(model, mvSaved, target, control) {
      gamma_node <- control$gamma_node
      random_walk <- nimble::nimbleFunctionList(nimble::sampler_BASE)
      random_walk[[1]] <- nimble::sampler_RW(model, mvSaved, target,
        list(scale = control$scale, adaptive = TRUE))
    },
    run = function() {
      if (model[[gamma_node]] == 0) {
        model$simulate(target)
        model$calculate(target)
        nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
      } else {
        random_walk[[1]]$run()
      }
    },
    methods = list(reset = function() { random_walk[[1]]$reset() }))
} else NULL

configure_omrf_mcmc <- function(model, E,
                                prior_type = "beta_binomial_independent",
                                proposal_scale = 0.20,
                                excluded_slab_gibbs = TRUE) {
  monitors <- c("interactions", "gamma", "rho")
  if (identical(prior_type, "beta_binomial_node")) {
    monitors <- c(monitors, "node_scale")
  }
  conf <- nimble::configureMCMC(model, monitors = monitors, enableWAIC = FALSE, print = FALSE)
  if (!prior_type %in% c("beta_binomial_independent", "beta_binomial_node")) {
    stop("Unknown OMRF NIMBLE prior type: ", prior_type, call. = FALSE)
  }

  # NIMBLE's generic configureRJ() did not update the coefficient when used
  # with these OMRF indicator nodes, and it explicitly does not support the
  # stochastic scales in the node-specific slab. Both targets therefore use an
  # exact fixed-dimensional slab_z augmentation. Keep the public proposal
  # scale control meaningful by applying it to the latent-slab RW samplers.
  for (node in paste0("slab_z[", seq_len(E), "]")) {
    conf$removeSamplers(node)
    edge <- sub("^slab_z\\[([0-9]+)\\]$", "\\1", node)
    if (excluded_slab_gibbs) {
      conf$addSampler(target = node, type = sampler_omrf_slab_z,
        control = list(scale = proposal_scale, gamma_node = paste0("gamma[", edge, "]")))
    } else {
      conf$addSampler(target = node, type = "RW",
        control = list(scale = proposal_scale, adaptive = TRUE))
    }
  }
  conf
}

omrf_nimble_defaults <- function(burnin = 2000L, iterations = 8000L,
                                 proposal_scale = 0.2) {
  list(niter = burnin + iterations, nburnin = burnin, thin = 1L, nchains = 1L,
       proposal_scale = proposal_scale)
}

run_omrf_nimble_rj <- function(stan_data, truth = NULL,
                               prior_type = "beta_binomial_independent",
                               burnin = 2000L, iterations = 8000L,
                               beta_binomial_a = 1, beta_binomial_b = 4,
                               start_label = NULL,
                               seed = 20260702L, run_fit = TRUE, compile = TRUE,
                               implementation = c("incident", "dense"),
                               excluded_slab_gibbs = TRUE, inits = NULL) {
  require_omrf_nimble()
  implementation <- match.arg(implementation)
  control <- omrf_nimble_defaults(burnin, iterations)
  control$implementation <- implementation
  control$excluded_slab_gibbs <- excluded_slab_gibbs
  control$sampler_type <- if (identical(prior_type, "beta_binomial_node")) {
    "latent_gaussian_scale_indicator"
  } else {
    "latent_cauchy_indicator"
  }
  if (!isTRUE(run_fit)) {
    return(list(
      status = "not_run",
      method = "nimble_rj",
      prior_type = prior_type,
      truth = truth,
      control = control,
      note = paste0(
        "Set RUN_FITS=TRUE to run the matching beta-binomial/",
        if (identical(prior_type, "beta_binomial_node"))
          "node-specific Gaussian" else "independent-Cauchy",
        " NIMBLE reference."
      )
    ))
  }

  setup_tick <- proc.time()[["elapsed"]]
  E <- stan_data$P * (stan_data$P - 1L) / 2L
  model <- omrf_nimble_model(
    stan_data,
    beta_binomial_a = beta_binomial_a,
    beta_binomial_b = beta_binomial_b,
    prior_type = prior_type,
    inits = inits %||% omrf_nimble_inits(stan_data, prior_type, seed,
      start_label = start_label),
    implementation = implementation
  )
  model$calculate()
  conf <- configure_omrf_mcmc(
    model,
    E,
    prior_type = prior_type,
    proposal_scale = control$proposal_scale,
    excluded_slab_gibbs = excluded_slab_gibbs
  )
  mcmc <- nimble::buildMCMC(conf)
  compiled <- NULL
  compile_seconds <- 0
  if (isTRUE(compile)) {
    compile_seconds <- system.time({
      compiled <- nimble::compileNimble(model, mcmc)
    })[["elapsed"]]
    run_mcmc <- compiled$mcmc
  } else {
    run_mcmc <- mcmc
  }
  setup_seconds <- proc.time()[["elapsed"]] - setup_tick

  stopifnot(control$thin == 1L, control$nchains == 1L)
  timing <- system.time({
    run <- run_mcmc_batches(run_mcmc, control$niter, control$nburnin,
                            as.integer(seed))
  })
  sample_matrix <- run$samples
  finalization_tick <- proc.time()[["elapsed"]]
  if (isTRUE(compile)) {
    # Release NIMBLE's compiled DLL and its finalizers before the worker is
    # recycled.  Leaving `compiled` alive until R teardown can produce the
    # harmless-looking but noisy `nimbleInternalFunctions` finalizer error.
    try(nimble::clearCompiled(compiled), silent = TRUE)
    rm(run_mcmc, compiled, mcmc, conf, model)
    invisible(gc(FALSE))
  }
  interaction_cols <- match_nimble_columns(sample_matrix, paste0("interactions[", seq_len(E), "]"))
  gamma_cols <- match_nimble_columns(sample_matrix, paste0("gamma[", seq_len(E), "]"))
  interaction_samples <- sample_matrix[, interaction_cols, drop = FALSE]
  gamma_samples <- sample_matrix[, gamma_cols, drop = FALSE]
  parameter_names <- paste0("interactions_0.", seq_len(E))
  colnames(interaction_samples) <- parameter_names
  colnames(gamma_samples) <- parameter_names
  ess <- coda::effectiveSize(interaction_samples)
  intervals <- t(apply(interaction_samples, 2, stats::quantile, probs = c(0.025, 0.975)))
  scale_samples <- NULL
  scale_summary <- NULL
  if (identical(prior_type, "beta_binomial_node")) {
    raw_names <- normalize_nimble_names(colnames(sample_matrix))
    scale_cols <- grep("^node_scale\\[[0-9]+\\]$", raw_names)
    if (length(scale_cols) != stan_data$P) stop(
      "Dependent OMRF NIMBLE output lacks the monitored node scales.",
      call. = FALSE)
    # PDMP's unconstrained scale coordinates are logarithms; retain both the
    # NIMBLE-native positive scales and their matched unconstrained values.
    scale_samples_raw <- sample_matrix[, scale_cols, drop = FALSE]
    scale_samples <- log(scale_samples_raw)
    scale_names <- sub("\\[", ".", sub("\\]", "",
      colnames(scale_samples_raw)))
    colnames(scale_samples) <- scale_names
    scale_ess <- coda::effectiveSize(scale_samples)
    scale_intervals <- t(apply(scale_samples, 2, stats::quantile,
      probs = c(0.025, 0.975)))
    scale_summary <- data.frame(parameter = scale_names,
      posterior_mean = colMeans(scale_samples),
      cri_low = scale_intervals[, 1L], cri_high = scale_intervals[, 2L],
      inclusion_probability = 1, ess = as.numeric(scale_ess),
      stringsAsFactors = FALSE)
  }
  interaction_summary <- data.frame(
    parameter = parameter_names,
    posterior_mean = colMeans(interaction_samples),
    cri_low = intervals[, 1L], cri_high = intervals[, 2L],
    inclusion_probability = colMeans(gamma_samples),
    ess = as.numeric(ess), stringsAsFactors = FALSE)
  posterior_summary <- if (is.null(scale_summary)) interaction_summary else
    rbind(interaction_summary, scale_summary)

  out <- list(
    status = "fit",
    method = "nimble_rj",
    prior_type = prior_type,
    start_label = start_label %||% "random_prior",
    truth = truth,
    control = control,
    fit = structure(list(
      stats = list(elapsed_time = unname(timing[["elapsed"]]),
        setup_time = unname(setup_seconds),
        compile_time = unname(compile_seconds),
        finalization_time = NA_real_,
        ct_ess = matrix(as.numeric(ess), nrow = 1L)),
      samples = interaction_samples,
      indicator_samples = gamma_samples,
      batches = run$batches,
      scale_samples = scale_samples,
      posterior_summary = posterior_summary
    ), class = "nimble_omrf_rj_result"),
    unc_names = posterior_summary$parameter,
    can_stick = c(rep(TRUE, E), rep(FALSE, nrow(posterior_summary) - E))
  )
  out$fit$stats$finalization_time <-
    proc.time()[["elapsed"]] - finalization_tick
  out$fit$stats$workflow_time <- setup_seconds +
    unname(timing[["elapsed"]]) + out$fit$stats$finalization_time
  out
}

normalize_nimble_names <- function(x) gsub(" ", "", x, fixed = TRUE)

match_nimble_columns <- function(samples, target_nodes) {
  matched <- match(normalize_nimble_names(target_nodes), normalize_nimble_names(colnames(samples)))
  if (anyNA(matched)) stop("Could not match NIMBLE sample columns: ", paste(target_nodes[is.na(matched)], collapse = ", "))
  matched
}

# ---- Starting values shared with the PDMP samplers ----------------------------

# NIMBLE initial values at the PDMP starting point: the posterior mode of the
# full model (every selectable coefficient included), from stan_mode()$constrained.
# RI-CLPM: the Stan and NIMBLE parameters have the same names; Stan's Cholesky
# factors of the correlation matrices are lower triangular and NIMBLE's
# dlkj_corr_cholesky factors upper triangular.
clpm_nimble_inits_from_mode <- function(mode) {
  m <- mode$constrained
  list(mu = as.numeric(m$mu), A = m$A, sd_ri = m$sd_ri, sd_w1 = m$sd_w1,
       sd_e = m$sd_e, Lcorr_ri = t(m$Lcorr_ri), Lcorr_w1 = t(m$Lcorr_w1),
       Lcorr_e = t(m$Lcorr_e))
}

# OMRF: thresholds and interactions from the mode of the independent-slab
# model, every edge included (gamma = 1, slab_z = interaction / slab sd), node
# scales 1, as for the PDMP samplers.
omrf_nimble_inits_from_mode <- function(mode, stan_data,
                                        prior_type = "beta_binomial_independent") {
  m <- mode$constrained
  P <- stan_data$P
  K <- stan_data$K
  thresholds <- matrix(0, P, K - 1L)
  start <- cumsum(c(1L, stan_data$seen - 1L))
  for (i in seq_len(P)) {
    n_i <- stan_data$seen[i] - 1L
    if (n_i > 0) thresholds[i, seq_len(n_i)] <- m$thresholds_0[start[i] - 1L + seq_len(n_i)]
  }
  node <- identical(prior_type, "beta_binomial_node")
  slab_sd <- if (node) stan_data$slab_base_scale else stan_data$prior_cauchy_scale
  inits <- list(rho = 0.2, gamma = rep(1L, length(m$interactions_0)),
                thresholds = thresholds,
                slab_z = as.numeric(m$interactions_0) / slab_sd)
  if (node) inits$node_scale <- rep(1, P)
  inits
}

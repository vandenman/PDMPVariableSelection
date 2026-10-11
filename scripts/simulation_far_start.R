# Time to the typical set from a start far from the posterior mode.
#
# Run from the project root, after scripts/simulation_omrf.R (whose mode-start
# fits give the reference posterior used by scripts/analyze_far_start.R):
#     Rscript scripts/simulation_far_start.R
#
# The OMRF data sets with p = 20, N = 200, random structure, and the
# independent prior are fitted again from the posterior mode plus `distance`
# in every unconstrained coordinate, in a random direction per coordinate (the
# same start for all methods). The PDMP samplers get no information from the
# mode: every coordinate starts with scale 1 and the Boomerang reference is
# centred at zero, and warmup adaptation proceeds as in the study. Only the
# warmup (PDMP) and burn-in (NIMBLE) are of interest, so the retained phases
# are short and both are recorded in fine batches. NIMBLE's reversible-jump
# sampler cannot start from distances of 15 or more (all log probabilities of
# its indicator updates are -Inf), hence distances 5 and 10.

# The OMRF study's definitions (data generation, fit_task()), without running
# it. run_tasks() marks its parallel workers with SIMULATION_WORKER before they
# source this file, so that they skip the run at the end.
is_worker <- identical(Sys.getenv("SIMULATION_WORKER"), "TRUE")
Sys.setenv(SIMULATION_WORKER = "TRUE")
source("scripts/simulation_omrf.R")
source("scripts/simulation_nimble.R")

# ---- Settings -----------------------------------------------------------------

distances <- c(5, 10)
workers <- 40
retained_T <- 10                     # PDMP: only the warmup (20 p = 400) matters
pdmp_settings$warmup_batches <- 1000L
nimble_burnin <- 10000
nimble_iterations <- 500
nimble_batches <- 2100L              # batches of 5 iterations
out_dir <- "simulation_results/far_start"

# ---- Design: the study's tasks (same data and seeds) at each distance ---------

far_tasks <- tasks[tasks$P == 20 & tasks$N == 200 & tasks$structure == "random" &
                     tasks$prior == "independent" &
                     (tasks$gradient %in% "full" | tasks$sampler == "nimble"), ]
far_tasks <- merge(far_tasks, data.frame(start_distance = distances), by = NULL)
far_tasks$id <- paste0(far_tasks$id, "_d", far_tasks$start_distance)

# ---- Run ----------------------------------------------------------------------

if (!is_worker) {
  run_tasks(far_tasks, setup_file = "scripts/simulation_far_start.R",
            out_dir = out_dir, workers = workers)
}

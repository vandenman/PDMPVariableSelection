# PDMPVariableSelection

Simulation code for the PDMP variable-selection experiments.

## Setup

Run all commands from this directory:

```sh
cd PDMPVariableSelection
```

On this machine, use the project-compatible R launcher (all examples below
assume the working directory shown above):

```sh
/home/don/bin/R -q
```

Restore the R package environment once before running simulations:

```sh
/home/don/bin/R -q -e 'renv::restore()'
```

The project uses the local development packages for this migration. Install the
R package into renv and select the local Julia package once:

```sh
/home/don/bin/R -q -e 'renv::install("/home/don/github/rpkgs/PDMPSamplersR", rebuild = TRUE)'
/home/don/bin/R -q -e 'PDMPSamplersR::check_for_julia_setup(); PDMPSamplersR::use_local_pdmpsamplers("/home/don/github/juliapkgs/PDMPSamplers.jl")'
```

Stan compilation, the persistent subset header, unconstrained parameter
mapping, and OMRF subsampling are all handled through exported
`PDMPSamplersR` interfaces. Content-addressed BridgeStan artifacts live in the
package cache rather than `stanmodels/`.

## Simulation entry points

The canonical entry points run the final ten-replication design by default.
Use `SIM_DRY_RUN=TRUE` to resolve and write a design manifest without fitting.

```sh
SIM_RETAIN_TRAJECTORY=TRUE Rscript scripts/simulation_clpm.R
SIM_RETAIN_TRAJECTORY=TRUE Rscript scripts/simulation_omrf.R
```

Dry-run outputs are written under the canonical namespaces:

- `simulation_results/clpm/final/`
- `simulation_results/omrf/final/`

## Run simulations

Fitting is the default. `SIM_WORKERS` is the worker-count setting
(`MIRAI_WORKERS` is retained only as a migration alias);
execution is dynamically scheduled and every condition is persisted atomically.

```sh
SIM_RETAIN_TRAJECTORY=TRUE Rscript scripts/simulation_clpm.R
SIM_RETAIN_TRAJECTORY=TRUE Rscript scripts/simulation_omrf.R
SIM_WORKERS=4 SIM_RETAIN_TRAJECTORY=TRUE Rscript scripts/simulation_clpm.R
SIM_WORKERS=4 SIM_RETAIN_TRAJECTORY=TRUE Rscript scripts/simulation_omrf.R
```

The OMRF manifest contains six PDMP routes for each prior: full and subsampled
Preconditioned Zig-Zag, full and subsampled Preconditioned BPS, and full and
subsampled Adaptive Boomerang. It also contains the NIMBLE reference rows at
`P = 10`. The accepted subsampled Adaptive Boomerang route uses
the existing person-node covariate-local harmonic certificate, `m = 1`, finite
anchor refresh distance `0.1`, no HCV, and no signed-harmonic research route.

Production PDMP condition files retain exact continuous-time posterior
summaries, inclusion probabilities, ESS, sampler counters, and—when the
command below enables retention—the materialized trajectory in a separate
atomic trace RDS. NIMBLE condition files preserve their sampled chains in the
ordinary result format and never require a PDMP trace.

OMRF production horizons are encoded by the declarative control table:
Preconditioned Zig-Zag uses warmup `5` and retained time `50`, Preconditioned
BPS uses warmup `5` and retained time `500`, and Adaptive Boomerang uses
warmup `50` and retained time `50`. Results are summary-only unless an
explicit research diagnostic option requests trajectory retention.

The only production command needed is therefore:

```sh
cd /home/don/hdd/surfdrive/Postdoc/Psychometrika_PDMP_Variable_Selection/PDMPVariableSelection
SIM_WORKERS=4 SIM_RETAIN_TRAJECTORY=TRUE Rscript scripts/simulation_clpm.R
SIM_WORKERS=4 SIM_RETAIN_TRAJECTORY=TRUE Rscript scripts/simulation_omrf.R
```

Use `SIM_METHODS` and the other small public filters only for a pilot. The
same command resumes compatible completed conditions and rejects results from
an incompatible manifest rather than overwriting them.

Results are saved incrementally in `simulation_results/<simulation>/final/conditions/`.
The top-level files below are small manifests containing the design and condition-file
paths; condition results are deliberately not combined in memory:

- `simulation_results/clpm/final/clpm_results.rds`
- `simulation_results/omrf/final/omrf_results.rds`

Each PDMP condition file contains summary-only continuous-time posterior
statistics and provenance. Full trajectories are research artifacts and are
not part of the production payload.

For short pilot runs, filter the grid with `SIM_METHODS`, `SIM_CONDITION_INDEX`, or `SIM_CONDITION_LIMIT`:

```sh
SIM_METHODS=PreconditionedZigZag SIM_CONDITION_LIMIT=1 \
SIM_WORKERS=1 Rscript scripts/simulation_clpm.R
```

## Methods

The CLPM simulation includes these full-gradient PDMP methods:

- `PreconditionedZigZag`
- `PreconditionedBPS`
- `AdaptiveBoomerang`

The CLPM crosses preconditioned Zig-Zag, preconditioned BPS, and Adaptive
Boomerang with
independent and dependent prior specifications. The dependent Adaptive
Boomerang route has been smoke-validated at `N = 100` and `N = 300` using the
beta-binomial model prior and a coherent exchangeable Gaussian slab over the
off-diagonal lag coefficients. Its full production rerun is recorded as a
follow-up pilot gate. The dependent specification uses sufficient statistics,
so subsampled CLPM methods are not part of the main simulation.

The CLPM analysis reports PIP AUC, threshold TPR/FPR, balanced accuracy,
completion/failure counts, and an explicit minimum-over-cross-lags ESS target
(`CLPM_FIXED_ESS_TARGET`, default `1000`). Current one-chain fixed-length
runs label R-hat and first ESS-target crossing time as unavailable rather than
silently imputing them; those two proposal commitments require a multi-chain
or ESS-stopping rerun.

The main OMRF PDMP comparison uses a beta-binomial model prior with an independent Cauchy slab and includes:

- `PreconditionedZigZag`
- `PreconditionedZigZag_subsampled`
- `PreconditionedBPS`
- `PreconditionedBPS_subsampled`
- `AdaptiveBoomerang`
- `AdaptiveBoomerang_subsampled`
The canonical OMRF design adds `nimble_rj` only at `P = 10`; it is scheduled
by the same resumable design but remains a distinct reference method.

All of these methods target the same posterior. The node-shared comparison
contains all six full/subsampled preconditioned routes, using the structured
log-linear Gaussian scale slab, plus its matching NIMBLE reference when
launched separately.
For node-shared subsampling, the analytic OMRF backend includes the
global-scale and node-scale Gaussian-prior gradients on their unconstrained
log coordinates; the likelihood residual envelope remains restricted to the
threshold and interaction coordinates.
The main OMRF design uses `P = 10, 20, 30` for PDMP and NIMBLE conditions.

The production PDMP method/prior matrix is explicit: each prior has the six
full/subsampled routes listed above. NIMBLE is physically separate. The
analysis writes `omrf_completion_summary.csv` with attempted,
completed, and failed counts by condition, method, and prior, including stored
failure messages.

The NIMBLE references use exact fixed-dimensional latent-slab augmentations:
a latent standard Cauchy for the independent Cauchy slab and a latent standard
normal multiplied by the stochastic edge scale for the node-shared Gaussian
slab. The persisted method identifier remains `nimble_rj` for compatibility;
`result$control$sampler_type` distinguishes `latent_cauchy_indicator` from
`latent_gaussian_scale_indicator`. `OMRF_NIMBLE_PROPOSAL_SCALE` sets the initial
adaptive random-walk scale for these latent slab variables.

The separate `simulation_omrf_bgms_appendix.R` experiment reproduces the `P = 10` main data-generating conditions and compares `zigzag_full` with `bgms_nuts` under a common independent Bernoulli/Cauchy prior. It is intentionally not part of the main factorial design.

Production subsampled OMRF methods use `m = 1`; the stochastic-gradient batch
size is recorded directly in the final design and is not inferred from a
nominal fraction.

## Script layout

Files directly under `scripts/` are the intended supplement-facing production
simulation, analysis, empirical-example, and shared-infrastructure entry
points. Temporary benchmarks, diagnostics, crossover screens, and failure
reproductions live under `scripts/Research/`; see its README for the boundary.
Superseded files remain under `scripts/deprecated/`.

## Useful Environment Variables

- `SIM_WORKERS`: number of isolated condition processes; default is `1`.
- `MIRAI_WORKERS`: temporary compatibility alias for `SIM_WORKERS`.
- `SIM_REPS`: number of replications; default is `10`.
- `SIM_SEED`: frozen base seed; default is `20260702`.
- `SIM_OUTPUT_DIR`: optional fresh output directory.
- `SIM_DRY_RUN`: resolve the design and skip fitting; default is `FALSE`.
- `SIM_RESUME`: reuse completed condition `.rds` files; default is `TRUE`.
- `SIM_PROGRESS`: show progress bars; default is `TRUE`.
- `SIM_METHODS`: comma-separated method filter.
- `SIM_COMPARISONS`: comma-separated comparison filter, for example `prior`.
- `SIM_PRIOR_TYPES`: comma-separated prior filter, for example `beta_binomial_node`.
- `SIM_N`, `SIM_P`, `SIM_K`: comma-separated integer filters for the corresponding design columns.
- `SIM_EDGE_STRUCTURES`: comma-separated OMRF edge-structure filter, for example `random,hub`.
- `SIM_REP_FILTER`: comma-separated replication-index filter.
- `SIM_MODE`: legacy `quick`/`full` compatibility switch; the default is `full`.
- `PDMP_CHAINS`: number of PDMP chains; default is `1`.
- `PDMP_SUBSAMPLE_SIZE`: generic override used by other scripts. The production
  OMRF script deliberately fixes `m = 1` for its validated subsampled methods.
- `PDMP_BETA_BINOMIAL_A` and `PDMP_BETA_BINOMIAL_B`: beta-binomial hyperparameters for the dependent CLPM model prior; defaults are 1 and 4.
- `PDMP_SLAB_SD` and `PDMP_SLAB_CORRELATION`: marginal standard deviation and exchangeable correlation of the dependent CLPM slab; defaults are 0.35 and 0.25.
- `OMRF_SLAB_BASE_SCALE`: fixed base scale multiplying the global and node scales in the dependent OMRF slab; default is 1.
- NIMBLE controls are kept in the model module and are not part of the
  production launcher interface.

## OMRF Analysis Plots

After running OMRF fits, regenerate summaries and plots with:

```sh
/home/don/bin/R -q -f scripts/analyze_simulation_omrf.R
/home/don/bin/R -q -f scripts/deprecated/analyze_omrf_bgms_appendix.R
```

Run the migration contract, compilation, gradient-closure, and envelope tests
with:

```sh
/home/don/bin/R -q -f tests/run_migration_tests.R
/home/don/bin/R -q -f tests/run_sampling_smoke_tests.R
/home/don/bin/R -q -f tests/run_omrf_nimble_regression_test.R
```

Main-study outputs are written to `simulation_results/omrf/final/plots/`.

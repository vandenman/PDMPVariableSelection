# PDMPVariableSelection

Code for the simulation studies, the empirical example and the figures of
*Accelerating Bayesian Variable Selection Using Piecewise Deterministic Markov
Processes*. The samplers are implemented in the Julia package
[PDMPSamplers.jl](https://github.com/vandenman/PDMPSamplers.jl) and used from R
through [PDMPSamplersR](https://github.com/vandenman/PDMPSamplersR).

## Requirements and setup

- R 4.6 and Julia 1.13.
- A C++ toolchain, for Stan (through BridgeStan) and NIMBLE.

Run everything from this directory. The R packages, including PDMPSamplersR,
are pinned in `renv.lock`:

```sh
Rscript -e 'renv::restore()'
```

PDMPSamplersR sets up its Julia packages (PDMPSamplers.jl and BridgeStan)
itself the first time it starts Julia. By default it uses a Julia project in
the user cache with the latest versions of these packages.

To use the versions with which the results were computed instead, point
PDMPSamplersR to the Julia project in `julia/` (`Project.toml` and
`Manifest.toml`); the pinned versions are then installed on first use. The
simulation studies run their fits in separate R processes, so set this in the
project's `.Rprofile`, which every process reads:

```r
options(pdmpsamplersr_config_julia_project_dir = normalizePath("julia"))
```

For a single interactive session,
`PDMPSamplersR::pdmpsamplers_setup(project_dir = "julia")` does the same.

## Scripts

All scripts are in `scripts/` and are run from the project root, e.g.
`Rscript scripts/simulation_omrf.R`. They can also be run step by step in an
interactive session. Each script documents its settings at the top.

| Step | Script | Produces | Manuscript |
| --- | --- | --- | --- |
| 1 | `plot_dynamics.R` | `dynamics.pdf` | Figure 1 |
| 2 | `simulation_clpm.R` | `simulation_results/clpm_study/` | RI-CLPM study |
| 3 | `analyze_clpm_study.R` | summaries, `clpm_ess_per_second_by_sampler.pdf` | Figure 2, RI-CLPM results |
| 4 | `simulation_omrf.R` | `simulation_results/omrf_study/` | OMRF study |
| 5 | `analyze_omrf_study.R` | summaries, `omrf_ess_per_second_by_sampler.pdf` | Figure 3, OMRF results |
| 6 | `analyze_convergence.R omrf` (and `clpm`) | time to convergence, `omrf_time_to_convergence.pdf` | Figure 4 |
| 7 | `simulation_far_start.R` (after step 4) | `simulation_results/far_start/` | Appendix F |
| 8 | `analyze_far_start.R` | `far_start_typical_set.pdf` | Figure F1 |
| 9 | `prepare_liss_health.R`, then `empirical_clpm.R` | `results/empirical_liss_health/` | Empirical example |
| 10 | `plot_empirical_clpm.R`, `table_empirical_clpm.R` | figure and appendix tables | Figure 5, Appendix E tables |

Shared code: `simulation_common.R` (fitting, starting values, online
summaries, parallel task runner), `simulation_nimble.R` (the NIMBLE
reversible-jump reference samplers), `figure_theme.R` (figure settings).
The Stan models are in `stanmodels/`.

### Computation

The studies run their fits in parallel R processes with
[mirai](https://mirai.r-lib.org/). The number of processes is set at the top
of each script (`workers`); the simulation studies and the empirical example
read it from the environment variable `SIM_WORKERS`. Every fit is saved as
soon as it finishes, so an interrupted study resumes where it stopped. The
computation time summed over all fits:

| Study | Fits | Default workers | Computation time |
| --- | ---: | ---: | --- |
| RI-CLPM | 1,440 | 10 | 76 hours (longest fit 26 minutes) |
| OMRF | 1,680 | 30 | 306 hours (longest fit 2.5 hours; NIMBLE 132 hours) |
| Distant start | 80 | 40 | under 20 minutes of wall time with 40 workers |
| Empirical example | 12 | 12 | 7 hours (longest fit 1.5 hours) |

The NIMBLE fits of the OMRF with p = 30 need about 7.5 GB of memory each;
choose `SIM_WORKERS` accordingly.

### Data of the empirical example

The empirical example uses the
[Health study of the LISS Core Study](https://www.dataarchive.lissdata.nl/study-units/view/12)
(<https://doi.org/10.17026/dans-ze3-5uk9>), which is available from the LISS
data archive after registration. The data cannot be shared here. Download the
English SPSS files of waves 2 to 7 (2008 to 2013; `ch08b` to `ch13g`) into a
folder `dataset/LISS/` next to this project, and run
`Rscript scripts/prepare_liss_health.R`; it writes
`../dataset/liss_health_panel.csv`, which `scripts/empirical_clpm.R` reads.

## Software versions and the reported results

The fits reported in the paper were computed with PDMPSamplers.jl and
PDMPSamplersR at the commits tagged `paper-rerun-2026-10` in both
repositories. The scripts here need the later versions pinned in `renv.lock`
and `julia/`, which simplify those packages without changing the samplers'
targets; see `software_versions/README.md` for the versions and the checks that were run.

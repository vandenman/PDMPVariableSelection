# Software versions

## Results reported in the paper (rerun of October 2026)

The simulation studies, the distant-start runs and the empirical example
were computed with

- PDMPSamplers.jl: commit 78254e3d46923169b5115c3627bb80a313955f0e,
  tag `psychometrika-rerun-2026-10`
- PDMPSamplersR: commit 57d32fe10da20fb39f184a24cbb96309d36d8bfe,
  tag `psychometrika-rerun-2026-10`
- Julia 1.13.0, R 4.6.1 (2026-06-24) "Happy Hop"
- Julia environment of PDMPSamplersR: `julia_Project.toml` and
  `julia_Manifest.toml` in this folder (PDMPSamplers.jl from the commit above;
  all other Julia packages as in the current `julia/` project)
- R packages: `renv.lock` of the paper repository's first commit (`bea6b5e`),
  with PDMPSamplersR from the commit above

## Current scripts

After the rerun, both packages were simplified (PDMPSamplers.jl PR #27,
PDMPSamplersR PR #14); obsolete options were removed, and the scripts in
`scripts/` were updated to the new interfaces. The scripts need

- PDMPSamplers.jl: commit 6591bba5b714f38b896bb7b9d612a388dde533f3
  (`main` after PR #27, tag `psychometrika-submission`), pinned in
  `../julia/Manifest.toml`
- PDMPSamplersR: commit 46d288de5ed450ab57534d9947f7cd2356019567
  (`main` after PR #14, tag `psychometrika-submission`), pinned in
  `../renv.lock`
- Julia 1.13, R 4.6.1

The simplifications were checked against fits of the rerun code: the first
round on PDMPSamplers.jl changed the random-number streams of some samplers,
and its fits agreed statistically (equal run times, five times longer than
in the study) with the rerun code; all later changes reproduced the fits of
the round before bit for bit (21 representative fits of both studies and the
distant-start runs, and the sampler paths of Figure 1). Rerunning the studies
with the current scripts therefore gives results that agree with the paper's
up to Monte Carlo error, not identical numbers.

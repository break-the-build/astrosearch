# Astrosearch research review

This is an evidence-backed review of the repository at commit
`fc0e95e99faddb92f20d5b94d6319ce4df33c9fc` (2026-10-10). The eight board
reports are sequential simulated reviews from one agent. They are not reviews
by the named scientists, nor independent human or multi-agent opinions.

## Current status

- Frozen baseline inventory and environment recorded.
- The known-answer release contains 200 calibration hosts and 200 blind hosts,
  each with a confirmed-planet curve and an inverted-flux control.
- All 600 required Rust searches completed. Input and output hashes match the
  manifests and sidecars.
- At a threshold selected on calibration controls, the production jackknife
  statistic recovered 87/200 blind known planets (43.5%, Wilson 95% CI
  36.8–50.4%) and triggered on 14/200 inverted controls (7.0%, CI 4.2–11.4%).
- The period grid hit its 250,000-trial cap on 583/600 searches (97.2%); median
  uncapped demand among those cases was 4.28 million trials. This is the
  strongest measured limitation.
- Three implementation cycles are documented: evaluation integrity, candidate
  lineage, and reproducible comparisons/telemetry. The third cycle adds new
  tools but awaits a full TLS cohort run and external hard-negative data. Ten
  cycles are not claimed; the inspected blind set is retired for tuning.
- Overall assessment: **Functional research prototype**. There is not enough
  evidence to claim survey completeness, calibrated astrophysical reliability,
  publication readiness, or a new discovery.

## Files

- [Baseline inventory](baseline/repository-inventory.md)
- [Environment and repository state](baseline/environment.md)
- [Baseline and post-change checks](baseline/baseline-tests.md)
- [Data release and integrity](baseline/data-manifest.md)
- [Eight simulated board reviews](boards/)
- [Adjudication and unresolved findings](adjudication/)
- [Prioritized improvement backlog](improvements/prioritized-backlog.md)
- [Completed cycle 1](improvements/cycle-01.md)
- [Completed cycle 2](improvements/cycle-02.md)
- [Cycle 3 implementation](improvements/cycle-03.md)
- [Recommendation status](improvements/implementation-status.md)
- [Hard-negative cohort protocol](improvements/hard-negative-protocol.md)
- [Final assessment and gates](final/)
- [Benchmark result summary](experiments/benchmark-results.md)

The frozen light curves and full per-case outputs are in the ignored local
directory `data/benchmarks/generated-20261009/` (about 2.4 GB). The committed
code contains the builder, runner, evaluator, and output auditor. Copy the
release directory to durable storage before cleaning ignored files.

## Evidence labels

**VERIFIED** means directly checked in repository files, manifests, hashes, or
executed results. **SUPPORTED** means consistent with authoritative documents
or papers, but not independently validated for every target here. **HYPOTHESIS**
means a testable explanation or proposed change. **REFUTED** means contradicted
by evidence. **UNKNOWN** means this review did not establish the claim.

## Scientific sources

- NASA [SPOC pipeline overview](https://heasarc.gsfc.nasa.gov/docs/tess/data_pipeline.html)
  and [TESS data documentation](https://heasarc.gsfc.nasa.gov/docs/tess/documentation.html).
- MIT's [Quick-Look Pipeline](https://tess.mit.edu/qlp-2/) and the [MAST QLP
  product page](https://archive.stsci.edu/hlsp/qlp), including its Sector 94
  change to TESS-Gaia Light Curve methodology.
- Hippke & Heller, [Transit Least Squares](https://arxiv.org/abs/1901.02015).
- Hippke et al., [Wōtan detrending benchmarks](https://arxiv.org/abs/1906.00966).
- Thompson et al., [Kepler DR25 completeness and reliability](https://doi.org/10.3847/1538-4365/aa984f).
- Giacalone et al., [TRICERATOPS validation of TESS candidates](https://arxiv.org/abs/2002.00691)
  and the updated [TRICERATOPS+ study](https://arxiv.org/abs/2508.02782).

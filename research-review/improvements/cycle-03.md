# Improvement cycle 03 — Make evaluations comparable and runs auditable

## Baseline and hypothesis

The review identified missing per-target runtime/build records, incomplete
stratification of the QLP product-era change, and no repeatable run command for
the existing TLS path across the frozen known-answer cohort. It also found
that synthetic vetting checks cannot estimate rates on real astrophysical
false positives.

**Hypothesis:** standardizing run metadata and evaluation inputs makes the next
real-data experiments reproducible; it does not, by itself, establish
scientific performance.

## Changes

- Added `scripts/run_metrics.py` to capture wall time, sampled process CPU,
  sampled peak RSS, source-tree fingerprint, source revision, executable hash,
  compiler versions, platform, Python version, and relevant package versions.
- Instrumented future known-answer and survey searches. Survey searches now
  write metrics to target manifests and automatically ingest successful
  detections into the shared candidate registry.
- Added `scripts/run_tls_benchmark.py`, which runs iterative TLS on the same
  frozen case light curves, calibrates the SDE trigger threshold using
  calibration controls only, and optionally compares paired blind outcomes
  against an astrosearch report.
- Added a hard-negative manifest freezer. The evaluator now reports
  known-false-positive trigger rates by externally supplied class, plus
  top-1/top-5/top-10/any-rank known-planet recovery.
- Benchmark manifests now record QLP extraction family and product sectors.
  The evaluator can backfill these fields from `case-metadata.csv` for the
  existing release and stratifies by extraction era, cadence, baseline,
  sector count, period/depth group, and truth class.
- Documented the hard-negative label protocol and TLS comparison command.

## Acceptance and remaining evidence

- Implementation is present; the historical 600-case release cannot gain
  accurate per-case runtime data retroactively. Its old executable hash also
  still has no archived source-commit mapping.
- The new TLS runner has not been run over all 600 cases. The current blind set
  is retired, so any comparison is retrospective and cannot tune production.
- No real hard-negative labels were added. The project needs externally
  dispositioned curves and references before reporting EB/blend/artifact rates.
- The existing search still caps a required 4.28M median period grid at 250k
  trials. The review has quantified the issue but has not replaced the search
  algorithm; no coarse-to-fine method is accepted without fresh injection and
  dense-oracle recovery evidence.
- GPU/Metal and empirical physical-vetting calibration remain deferred pending
  a new, appropriately labeled cohort and a search algorithm with verified
  resolution.

The next work item is to run the new comparison and stratification tools on
the frozen release for descriptive analysis, then collect new host-level
calibration/blind samples before any tuning.

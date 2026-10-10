# Improvement cycle 01 — Make benchmark evaluation match the deployed tier

## Baseline and hypothesis

The generic benchmark evaluator calibrated raw `score`, while the current
candidate tier and `data/calibration/detection_threshold.json` use
`jackknife_score`. The separate calibrator silently skipped cases without
outputs and selected the first matching candidate rather than the highest
matching statistic.

**Hypothesis:** calibrating/evaluating the production jackknife statistic,
requiring all cases used by the protocol, and taking the maximum among matching
candidate ephemerides will produce a reproducible measure of the decision rule
actually used by the project.

## Changes

- `scripts/evaluate_benchmark.py`: added a statistic selector; the CLI defaults
  to the production `jackknife_score`; threshold selection, trigger counts, and
  ephemeris recovery now use that statistic. The generic helper retains raw
  score as its backward-compatible default.
- `scripts/calibrate_detection_threshold.py`: only calibration controls and
  blind cases are required; omitted calibration positives are explicit. Missing
  required outputs now stop evaluation instead of shrinking a denominator.
  Planet recovery now uses the maximum statistic among all matching candidates.
- Added `scripts/audit_known_answer_benchmark.py` for coverage/provenance/hash
  verification and period-grid cap reporting.
- Added `scripts/compare_detection_statistics.py` for paired host bootstrap and
  exact McNemar comparisons of frozen evaluation reports.
- Added regression tests for jackknife calibration, matching-candidate choice,
  missing required outputs, and paired outcomes.

## Acceptance and result

- Required manifest: 600/600 exact cases; all inputs and candidate outputs
  match their hashes; one binary and search configuration across the run.
- Baseline Rust tests: 62 passed. Python unit tests: 57 passed before changes;
  61 passed after this cycle's evaluator changes (before the separate registry
  cycle).
- The calibrated production statistic threshold is 22.4876456 from 10/200
  calibration-control triggers.
- Blind jackknife result: 87/200 recovery (43.5%, Wilson 95% CI 36.8–50.4%);
  14/200 inverted-control triggers (7.0%, CI 4.2–11.4%).
- Against the separately calibration-thresholded raw score, jackknife adds
  26 recoveries and loses none among positives; paired difference +13.0 points
  (95% bootstrap +8.5 to +18.0, McNemar p=2.98e-8). Inverted triggers rise
  from 8 to 14; paired difference +3.0 points (95% bootstrap +0.5 to +6.0,
  McNemar p=0.070).

This cycle changes evaluation integrity, not the Rust search output. It exposes
the measured trade-off instead of changing the threshold to favor a desired
result. The known-answer blind split is now retired from further tuning.

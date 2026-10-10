# Prioritized research and engineering backlog

| Priority | Item | Why | Acceptance evidence |
| --- | --- | --- | --- |
| P0 | Resolve or quantify period-grid cap | **Quantified; unresolved.** 583/600 cases needed more than the 250k cap; median demand 4.28M | Still requires a new locked injection/host holdout comparing a coarse-to-fine/event search against dense local oracles. Current inspected blind set is retired. |
| P0 | Add hard astrophysical negatives | **Evaluation path implemented; data missing.** `known_false_positive` is reported by subtype and the curated-manifest builder validates label references/checksums | Add externally dispositioned EB, blend, variable, and artifact cases; report per-class trigger/vetting rates with CIs. |
| P1 | Compare Rust box, TLS, SPOC/QLP | **Runnable paths implemented; comparative run pending.** TLS benchmark runner added; `tess_survey.py` supports QLP/SPOC and TLS | Run on a common frozen host set with thresholds calibrated only on calibration controls; retain pipeline/product versions and paired outcomes. |
| P1 | Archive exact release and build identity | **Future-run identity implemented; historic release incomplete.** Code, compiler, executable, platform and dependency fingerprints are now recorded | Recreate the old binary-to-source mapping or freeze a new benchmark release, then archive the data and build under a DOI. |
| P1 | Record per-target operational metrics | **Implemented for new runs.** Search runner records wall time, sampled CPU/RSS, configuration, build identity, warnings/log path and exit status | Re-run on a new dataset to populate metrics; old benchmark outputs remain without runtime telemetry. |
| P1 | Stratify QLP processing-family changes | **Manifest/evaluator support implemented.** Builders and evaluator can separate legacy QLP from the S94+ TGLC-based family, cadence, baseline, sector count and target strata | Re-score the frozen release from case metadata, then confirm with product-specific sensitivity runs. |
| P1 | Expand injection/recovery surface | **Mechanism exists; adequate survey surface not measured.** Real-curve harness supports randomized depth, period, phase, duration and impact parameter, finite exposure, gaps and season holdout | Run a preregistered multi-target injection grid and publish recovery contours with host-clustered intervals. |
| P2 | Improve candidate ranking and depth uncertainty | **Top-k recovery added; uncertainty calibration pending.** Evaluator now reports top-1/5/10/any-rank | Validate depth uncertainty and sector consistency against a fresh labeled set; never use the current blind split for tuning. |
| P2 | Validate physical vetting modules on real labels | **Pending external labels.** Existing tests cover synthetic scenarios | Run blinded vetting on archived dispositions and assess false-positive rates, uncertainty coverage, and reviewer agreement. |
| P3 | Evaluate GPU/Metal | **Deferred.** Period-grid sampling dominates correctness risk | After the search redesign, compare GPU/Metal with a 64-bit CPU oracle for fidelity, throughput, energy and memory. |

Do not tune the next release against the now-inspected blind cases. Keep the
current threshold and search settings fixed for reporting; use new targets or
registered injections for the next experiment.

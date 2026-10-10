# Board 8 — Publication, reproducibility, and open science

**Sequential simulated review; not an independent journal review.**

## Finding

**Rating: 2/5 — reproducible local evaluation is emerging; publication package
and novel result are not established.**

Strengths include source code, tests, frozen NASA catalogue provenance,
case-level hashes, a label-free run manifest, calibration/blind splits, and a
machine-readable evaluator. The generated release is currently ignored by
Git and occupies about 2.4 GB. The search executable is not archived, and the
sidecar has a binary hash but no source commit/build identity. The run did not
record per-case elapsed time. An independent second person has not repeated
the benchmark.

The science question is presently detector behavior on published positives and
inverted controls. It does not establish a new exoplanet, survey completeness,
astrophysical reliability, or a new method result against a strong baseline.
The known-answer sample has now been inspected; it cannot serve as an untouched
blind test for future tuning.

## Publication prerequisites

- Archive input release, truth/run manifests, build source and exact executable
  under an immutable, citable DOI with usage/license and checksums.
- Document source product versions, processing-family changes, selection
  filters, all exclusions, and the distinction between calibration and blind
  metrics.
- Add hard-negative and injection/recovery cohorts; compare with a defensible
  independent baseline at matched false-alarm rates.
- Report selection effects, confidence intervals, period/depth/cadence strata,
  period-grid cap rate, and limitations.
- Obtain external reproduction and candidate-specific independent validation
  before describing any object as a discovery.

**Gate:** publication readiness is **INSUFFICIENT EVIDENCE**. A methods or
software note might become appropriate after the benchmark is archived and
the search-space cap is addressed; a discovery paper is not supported by this
review.

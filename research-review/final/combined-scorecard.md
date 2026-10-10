# Combined scorecard

Scores are the primary agent's evidence-based simulated ratings, not actual
opinions from named individuals or independent agents. The weighted result is
descriptive only; it cannot override release gates.

| Board | Weight | Rating / 5 | Evidence summary |
| --- | ---: | ---: | --- |
| Scientific review | 20% | 3.0 | Plausible light-curve screen; not a validated discovery pipeline |
| Software engineering | 15% | 3.5 | Modular Rust/Python code and strong tests; runtime/build provenance gaps and period cap |
| Statistics and validation | 20% | 2.5 | Real 200-host holdout, but 43.5% recovery and only inverted nulls |
| Data integrity | 15% | 3.5 | Strong checksums/provenance; QLP method changes and no pixel-level benchmark |
| Algorithm research | 10% | 2.5 | Efficient shared preprocessing; 97.2% trial cap and no matched TLS/SPOC evaluation |
| Adversarial red team | 10% | 2.0 | No hard-negative cohort; several astrophysical alternatives remain untested |
| Observational confirmation | 5% | 2.0 | Vetting code exists; no candidate in this review received external confirmation |
| Publication/reproducibility | 5% | 2.0 | Local hashes and scripts; release data/executable not archived and no independent rerun |

**Weighted simulated score: 2.8/5.** The gate result is more important than
this aggregate: no scientific-validation or external-research gate passes.

## Critical or high-impact findings

1. The period search undersamples its declared timing resolution on 583 of 600
   targets; this is a verified limitation, not a hypothetical optimization.
2. False-alarm calibration uses inverted light curves only. It is not the
   astrophysical false-positive rate.
3. Recovery is conditional on a discovery-selected confirmed-planet sample and
   is not a survey completeness estimate.
4. Detection-tier evaluation was inconsistent across scripts; cycle 1 fixed
   the generic evaluator, missing-case handling, and matched-candidate
   aggregation.
5. The benchmark products are not yet an immutable public reproducibility
   package; binary hash is known but the executable/source mapping is not.

## Overall verdict

**Functional research prototype.** The project is useful for candidate triage
and research experiments. Claims of calibrated survey sensitivity, reliable
new-planet discovery, and publication readiness require the P0/P1 backlog.

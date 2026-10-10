# Evidence adjudication

## Decision

**Overall verdict: Functional research prototype.** The software can process
TESS light curves and recover a substantial fraction of an easy,
discovery-selected confirmed-planet sample. The measured sensitivity is not yet
high enough, the grid is capped on nearly every target, and false-alarm
controls are too narrow to support a discovery claim or survey reliability
claim.

No simulated reviewer vote can override those limits. No board reported a
verified fundamental scientific error in the basic weighted box search, but
the current search resolution and validation coverage block research-release
claims beyond prototype-level candidate triage.

## Best-supported positive case

- Search methods and output are explicit; the pipeline does not claim that a
  high score confirms a planet.
- Rust tests cover numerical, cadence, mask, alias, and serial/parallel
  invariants; Python tests cover catalog, FAP, vetting, and resume logic.
- A 200-host calibration split and a disjoint 200-host blind split were
  actually searched. The run has exact case-level inputs/outputs and matching
  hashes.
- The jackknife statistic has a reproducible paired advantage over raw score
  for recovery on this blind sample.

## Strongest case against a scientific-acceptance claim

- Only 87/200 blind confirmed planets pass the calibrated production tier.
- Blind inverted-control trigger rate is 14/200, with a wide interval.
- No known false-positive/eclipsing-binary cohort was evaluated.
- 583/600 searches hit the 250,000-period cap, whose warnings explicitly state
  narrow peaks may be missed.
- The benchmark is selected from already-confirmed systems; it cannot estimate
  population completeness.
- No new signal underwent independent pipeline, pixel-source, or observing
  confirmation in this review.

## Reconciled metric disagreement

Raw-score evaluation and jackknife-tier evaluation measure different rules.
Each threshold was selected on the same 200 calibration controls at a 5%
target. On the blind sample, jackknife recovered 87/200 and raw score 61/200.
The paired recovery difference was +13 points (95% host-bootstrap interval
+8.5 to +18.0; exact McNemar p=2.98e-8). The corresponding inverted-control
trigger rates were 7% and 4% (difference +3 points, interval +0.5 to +6.0,
McNemar p=0.070). The evidence supports using the jackknife statistic for
candidate tiering on this cohort, with a possible false-alarm cost. It does
not prove that either threshold is universally calibrated.

## Independent-review limitation

The eight reports are sequential simulations by one agent. A separate review
context/agent was not available under the active tool policy. Therefore there
was no independent initial review, no inter-board blind period, and no real
cross-examination. Their scores are an organized internal assessment, not
independent expert consensus.

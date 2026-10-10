# Conflicting findings

## Detection statistic

Raw score has higher held-out specificity in this sample (8/200 inverted
triggers) but lower recovery (61/200). The jackknife statistic discounts a
candidate dominated by one event; it recovers 87/200, with 14/200 inverted
triggers. Paired data support a recovery gain; the paired increase in false
alarms is suggestive but not conventionally significant at 5% (McNemar
p=0.070). Do not summarize this as a free improvement: the operating point
changes.

## Current threshold versus blind estimate

`data/calibration/detection_threshold.json` stores the all-control deployment
threshold 23.3025. That threshold uses all 400 inverted controls after the
benchmark analysis. The independent blind-performance estimate instead uses
the calibration-only threshold 22.4876 and must remain the reported holdout
result. The deployment threshold is not a second blind-test result.

## Grid coverage versus measured recall

The benchmark gives an empirical recovery estimate for the capped search that
was run. The high cap frequency is a plausible contributor to misses, but
results for capped and uncapped cases cannot establish causation: only eight
blind targets are uncapped, and target characteristics differ. A full-grid or
coarse-to-fine paired rerun is required to attribute the misses.

## Inverted controls versus real false positives

Inversion is useful for instrumental/noise controls and is supported by the
Kepler DR25 reliability literature. It does not simulate the full astrophysical
false-positive population. These results are specifically the inverted-control
trigger rate, not an all-cause FPR.

## Scientific support versus publication readiness

The algorithm is scientifically plausible as a screening method, but current
data do not support publication-ready completeness, reliability, or a discovery
claim. Those are different acceptance questions.

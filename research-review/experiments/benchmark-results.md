# Known-answer benchmark results

## Protocol

Release: `data/benchmarks/generated-20261009/`. The 200 calibration hosts and
200 blind hosts are disjoint. Each host has one confirmed-planet QLP curve and
one matched inverted-flux control. The detection threshold was selected only
from calibration inverted controls at a 5% target; the blind set did not select
the threshold. The production statistic is the event jackknife score.

The blind set has been examined and is no longer an untouched holdout for
future tuning. Keep these numbers as the result for this frozen release and
reserve a new sample or preregistered injections for further changes.

## Primary operating point

| Metric | Result | 95% Wilson interval |
| --- | ---: | ---: |
| Calibration inverted-control trigger | 10/200 = 5.0% | 2.7–9.0% |
| Blind known-planet recovery | 87/200 = 43.5% | 36.8–50.4% |
| Blind inverted-control trigger | 14/200 = 7.0% | 4.2–11.4% |

The calibrated jackknife threshold was 22.4876. The paired benchmark contains
200 independent blind planet hosts and 200 corresponding controls. The
worst-case Wilson half-width at n=200 is about 6.9 percentage points; the
observed recovery interval half-width is 6.8 points.

At any candidate rank/score, 160/200 known planets had a matching periodic
candidate. Only 87/200 passed the jackknife tier. This points to losses both in
search recovery and in candidate ranking/tiering.

## Candidate ranking and QLP extraction-era strata

On the same frozen blind positives, the top-ranked transit candidate matched
the catalog ephemeris for 140/200 hosts (70.0%; 95% Wilson CI 63.3–75.9%). A
match occurred in the top five for 158/200 (79.0%; 72.8–84.1%) and top ten for
160/200 (80.0%; 73.9–85.0%). The total score threshold still admitted 87/200;
ranking and calibrated tier selection measure different things.

| QLP product-era group | Blind planet recovery | Blind inverted-control trigger |
| --- | ---: | ---: |
| Legacy QLP only | 59/140 = 42.1% (34.3–50.4%) | 8/140 = 5.7% (2.9–10.9%) |
| Mixed legacy and S94+ products | 28/60 = 46.7% (34.6–59.1%) | 6/60 = 10.0% (4.7–20.1%) |

There were no blind hosts with only the S94+ extraction family, so this does
not estimate performance on S94+ products alone. The groups are small and
selected, the products within a stitched light curve can span multiple
sectors, and the differences are descriptive rather than causal. The
machine-readable re-evaluation is
`data/benchmarks/generated-20261009/evaluation-jackknife-by-product.json`.

## Depth strata

| Approximate depth | Recovered | Rate | Wilson 95% interval |
| --- | ---: | ---: | ---: |
| Deep | 46/61 | 75.4% | 63.3–84.5% |
| Medium | 36/113 | 31.9% | 24.0–40.9% |
| Shallow | 4/26 | 15.4% | 6.2–33.5% |

The sample size is small in each stratum and the positives are selected from
previously confirmed objects. Treat these as conditional benchmark estimates,
not completeness contours.

## Paired statistic comparison

At thresholds independently chosen from the calibration controls (both 5%
training target), jackknife versus raw score gave:

| Blind measure | Jackknife | Raw score | Paired difference (jackknife − raw) |
| --- | ---: | ---: | ---: |
| Planet recovery | 87/200 = 43.5% | 61/200 = 30.5% | +13.0 pp (95% bootstrap +8.5 to +18.0; McNemar p=2.98e-8) |
| Inverted-control trigger | 14/200 = 7.0% | 8/200 = 4.0% | +3.0 pp (95% bootstrap +0.5 to +6.0; McNemar p=0.070) |

This supports the jackknife score as a better recovery ranker on this release,
with a measured increase in inverted-control triggers. It does not establish
all-cause false-positive probability.

## Search resolution and provenance

- 600/600 required searches completed; no failed/missing cases.
- Input and output hashes matched for all 600 cases.
- The recorded binary hash is `638f51dece3099b37bafed67ea88ea60753137d6f10fcb8b6e04666a16760c84`.
- 583/600 searches (97.2%) hit the 250,000-trial cap. Capped searches requested
  599,204–6,558,648 trials; median requested count was 4,280,132.
- 17 searches were uncapped. Only eight blind cases were uncapped, so a
  capped-versus-uncapped recovery comparison is too small and confounded to
  attribute causal effects.
- Per-case runtime was not recorded. The output directory contains 183
  unreferenced sidecars from earlier runs; the exact 600-row manifest controls
  evaluation and excludes them.
- No hard astrophysical negatives were present; known-false-positive trigger
  rates remain unmeasured.

## Reproduction

```sh
.venv/bin/python scripts/audit_known_answer_benchmark.py \
  --run-manifest data/benchmarks/generated-20261009/run-manifest.csv \
  --truth-manifest data/benchmarks/generated-20261009/truth-manifest.csv \
  --provenance data/benchmarks/generated-20261009/provenance.json \
  --output data/benchmarks/generated-20261009/audit.json

.venv/bin/python scripts/evaluate_benchmark.py \
  --manifest data/benchmarks/generated-20261009/truth-manifest.csv \
  --output data/benchmarks/generated-20261009/evaluation-jackknife.json \
  --target-false-alarm-rate 0.05 --statistic jackknife_score

.venv/bin/python scripts/evaluate_benchmark.py \
  --manifest data/benchmarks/generated-20261009/truth-manifest.csv \
  --output data/benchmarks/generated-20261009/evaluation-score.json \
  --target-false-alarm-rate 0.05 --statistic score

.venv/bin/python scripts/compare_detection_statistics.py \
  data/benchmarks/generated-20261009/evaluation-jackknife.json \
  data/benchmarks/generated-20261009/evaluation-score.json \
  --output data/benchmarks/generated-20261009/statistic-comparison.json
```

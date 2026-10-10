# Board 3 — Discovery validation and statistics

**Sequential simulated review; not an independent external audit.**

## Finding

**Rating: 2.5/5 — useful first empirical calibration, insufficient for a
survey-level reliability claim.**

The release contains 200 host stars for calibration and 200 disjoint blind
hosts. Each host has one catalogued positive and one inverted control. Threshold
selection uses only the 200 calibration controls. The calibrated jackknife
threshold is 22.4876; 10/200 calibration controls trigger (5.0%, Wilson 95% CI
2.7–9.0%). On the blind split:

- 87/200 confirmed planets recover (43.5%, Wilson 95% CI 36.8–50.4%).
- 14/200 inverted controls trigger (7.0%, Wilson 95% CI 4.2–11.4%).
- 160/200 planets have a matching candidate at some score, before the detection
  threshold; score/ranking/tiering therefore loses additional recoveries.
- The 200-host sample has a worst-case Wilson half-width of about 6.9 points.

At the same 5% calibration target, jackknife thresholding recovered 26 blind
planets that raw-score thresholding missed and lost none of raw-score's 61
recoveries. The paired host bootstrap difference was +13.0 points (95%
percentile interval +8.5 to +18.0); exact McNemar p=2.98e-8. Jackknife triggered
on 14 inverted controls versus raw score's 8; paired difference +3.0 points
(95% interval +0.5 to +6.0), McNemar p=0.070. This favors the production
jackknife statistic for recovery, with a possible false-alarm cost; it does not
validate its probability calibration.

The null class is only inverted-flux curves. There is no known-false-positive
cohort in this release. A 7% held-out inverted-control rate is compatible with
the calibration target given its uncertainty, but says nothing about the FPR
on eclipsing binaries, blends, or astrophysical false positives.

## Recommendations

- Use host-level paired analyses for algorithm comparisons. Keep positive and
  inverted-control rows paired by host.
- Add independent real hard negatives and realistic injections over depth,
  duration, period, impact parameter, transit phase, cadence, and gaps.
- Measure recall at rank 1/top-k, at fixed calibration thresholds, and across
  strata. Report confidence intervals and sample sizes for every stratum.
- Retire this blind split now that it has been evaluated. Any future code or
  threshold tuning requires a new blind set or registered injection suite.
- Keep the meaning of score distinct from posterior probability, FAP, or
  planet-confirmation confidence.

## Status

**VERIFIED:** thresholds and Wilson intervals were computed from frozen
per-host outputs. **UNKNOWN:** the survey-wide completeness surface and the
empirical astrophysical false-positive rate.

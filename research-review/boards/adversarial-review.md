# Board 6 — Adversarial red team and falsification

**Sequential simulated review; not an independent external audit.**

## Strongest alternative explanations for a high-scoring transit

1. A one- or two-event instrumental excursion is being folded into a plausible
   period. The project's jackknife tier addresses single-event domination,
   but the transit score remains heuristic and period-grid capping can distort
   peak selection.
2. The signal is an eclipsing binary, a background/nearby eclipsing binary,
   or a blended source. Inverted curves do not simulate these cases.
3. Sector offsets or processing changes create an apparent periodic dip.
4. A known TOI/TCE/planet or a published candidate is being rediscovered.
5. A stitched product or time conversion is wrong at a sector boundary.

## Finding

**Rating: 2/5 — useful checks exist, but the sample does not falsify the main
astrophysical alternatives.**

**VERIFIED:** the benchmark preserves opaque run IDs and keeps its truth table
out of the run manifest; all expected inputs/outputs hash correctly. Existing
vetting code checks secondary eclipses, odd/even depths, sector consistency,
Gaia neighbors, centroid offsets, catalog status, and empirical FAP. These
tools were not all run against every known-answer case as a validation study.

**VERIFIED LIMITATION:** there are no catalogued eclipsing-binary/known-false-
positive rows in the evaluation; the known-FP rate is therefore not estimable.
583 searches also hit the period cap. **UNKNOWN:** whether the resulting
candidate list is robust to independent extraction, pixel-level localization,
or an independent human analysis.

## Required red-team experiments

- Preserve real eclipsing binaries and variable stars as labeled hard
  negatives; evaluate classification separately from signal detection.
- Use sector permutations/time shifts and injected transits to challenge
  sector-offset artifacts.
- Inspect highest-score false triggers on original pixel data and neighboring
  Gaia sources.
- Compare candidate ephemerides against TOI/CTOI/TCE and confirmed-planet
  catalogues at the catalogue epoch.
- Rerun a subset with a second detrending and independent TLS/SPOC search.
- For every claimed candidate, document the strongest non-planet explanation
  and the observation that would distinguish it.

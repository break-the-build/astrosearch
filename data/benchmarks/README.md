# Known-answer benchmark

This benchmark measures recovery of published, confirmed transiting planets in
real TESS QLP light curves, and measures false triggers on paired inverted-flux
controls made from the same light curves. It is designed for repeatable
algorithm comparisons and score-threshold calibration. It does **not** estimate
the occurrence rate of planets, survey-wide completeness, or planet-confirmation
reliability.

## Known objects

The sample is built from the NASA Exoplanet Archive `pscomppars` table, filtered
to confirmed transiting planets discovered by TESS, with a TIC identifier,
period 0.15–30 days, a transit midpoint and duration. The builder retains one
planet per TIC so target stars, rather than planets in multi-planet systems,
are the independent sampling units. Examples of confirmed TESS planets the
Archive recognizes include **Pi Mensae c**, **LHS 3844 b**, and **TOI-700 d**;
the generated frozen catalog snapshot records the actual objects selected for
each benchmark release.

NASA's TAP service documents the archive tables and supports reproducible CSV
queries. [TAP guide](https://exoplanetarchive.ipac.caltech.edu/docs/TAP/usingTAP.html),
[PS/PSCompPars column definitions](https://exoplanetarchive.ipac.caltech.edu/docs/API_PS_columns.html).
NASA's TESS releases describe the discoveries of [Pi Mensae c and LHS 3844 b](https://science.nasa.gov/universe/exoplanets/nasas-tess-rounds-up-its-first-planets-snares-far-flung-supernovae/)
and the [TOI-700 system including planet d](https://www.nasa.gov/universe/nasa-planet-hunter-finds-its-1st-earth-size-habitable-zone-world/).
Light curves are downloaded from MAST through Lightkurve's TESS QLP product
search. [Lightkurve/MAST data guide](https://lightkurve.github.io/lightkurve/tutorials/1-getting-started/searching-for-data-products.html).

## Statistical design

The default build requests 200 distinct host stars for threshold calibration
and another 200 distinct host stars for blind evaluation. Every host supplies
one positive known-planet light curve and one inverted-flux control. The split
is fixed by a documented seed and host TIC, before any search is run. A target
star appears in only one split. With 200 blind hosts per class, the worst-case
95% Wilson interval half-width for a binary rate is about **6.9 percentage
points**. Reports always show the observed interval; fewer cases are marked
underpowered. For ±5 percentage points at 95% confidence, increase each split
to at least 385 independent hosts (`--cases-per-split 385`).

Threshold calibration uses only the calibration split's inverted controls. It
chooses the lowest score threshold whose observed calibration trigger rate is
at or below the requested rate (5% by default). The threshold is then frozen
and evaluated once on blind positives and inverted controls. Do not change
search logic or thresholds after inspecting the blind results; create a new
benchmark release or reserve a new blind sample for further tuning.

The inverted controls preserve each target's cadence, gaps, uncertainties,
stellar variability, and instrumental background while reversing flux
excursions. They are useful empirical controls, following the use of inverted
light curves in Kepler DR25 reliability studies, but they do not represent the
full population of eclipsing binaries, blends, or other astrophysical false
positives. See [Thompson et al. 2018, Kepler DR25 catalog reliability and
completeness](https://doi.org/10.3847/1538-4365/aa984f).

## Build and run

The normal build requires network access to the NASA Exoplanet Archive TAP
service and MAST. It downloads, checks, and freezes the catalog CSV, light
curves, metadata, and SHA-256 hashes under an ignored generated-data directory.
Do not use a pilot build as a statistically powered calibration.

```sh
uv sync
uv run python scripts/build_known_answer_benchmark.py \
  --output data/benchmarks/generated \
  --cases-per-split 200 \
  --seed 20261009 \
  --workers 4 \
  --timeout 90
cargo build --release
uv run python scripts/run_known_answer_benchmark.py \
  --manifest data/benchmarks/generated/run-manifest.csv \
  --threads 1 --workers 24 --transit-search-passes 1
uv run python scripts/evaluate_benchmark.py \
  --manifest data/benchmarks/generated/truth-manifest.csv \
  --output data/benchmarks/generated/evaluation.json \
  --target-false-alarm-rate 0.05
```

`run-manifest.csv` has only opaque case IDs, mission/split labels, input paths,
and output paths. Keep `truth-manifest.csv`, `case-metadata.csv`, and
`nasa-pscomppars-tess-transiting.csv` out of version control until the blind
evaluation is retired. `case-metadata.csv` contains planet names, host IDs,
ephemerides, product IDs, and checksums; `provenance.json` records the exact
archive query, response checksum, release time, seed, Lightkurve version,
sample size, and known limitations. The evaluator reports recovery by period
and estimated-depth strata, with per-case cadence, depth proxy, stellar
temperature, and magnitude metadata. `run_known_answer_benchmark.py`
is resumable and reuses a result only when the input, executable, and search
configuration hashes match. `--rerun` forces a fresh search.

The benchmark has one known planet per host, so its recovery run uses one Rust
transit-search pass; this evaluates recovery of that labeled signal without
spending compute on secondary-planet searches. The normal survey default stays
multi-pass. Run targets in parallel with `--workers`; set `--threads 1` per
process to keep total CPU use bounded. The label-free run manifest contains
calibration controls and all blind cases; calibration positives remain in the
separate truth manifest and are not searched. The per-case sidecars record the
exact binary hash and search options.

For publications, cite the NASA Exoplanet Archive PS Composite Parameters
dataset DOI [10.26133/NEA13](https://doi.org/10.26133/NEA13) and the relevant
MAST product documentation; record the exact QLP product IDs in the release
manifest.

Use `--cases-per-split 20 --pilot` only to check setup and product access. A
pilot is explicitly marked `underpowered`; it must not be cited as detection
accuracy. The normal builder fails instead of silently shrinking the requested
sample. Downloads use four isolated workers by default; lower `--workers` if
MAST rate-limits requests. `--timeout` bounds network requests. Increase to
385 per split for approximately ±5 percentage-point worst-case precision. As
the NASA catalog and MAST products change over time,
preserve each released `generated/` directory outside the ignored workspace,
alongside source and executable hashes, to make future comparisons exact.

## Interpretation and future extensions

The known-planet set is discovery-selected: it is weighted toward planets that
were already easy enough to find and confirm. The measured positive rate is
therefore **conditional recovery on this benchmark**, not completeness for
TESS's full target population. The matched inverted controls are paired by
host, but the per-class Wilson intervals are descriptive binomial intervals;
for comparing two algorithm versions, use paired host-level tests or a
host-cluster bootstrap. Report TESS sectors, cadence, period, transit duration,
depth proxy, stellar type, and magnitude strata so pooled rates do not hide
weak regimes.

Keep known eclipsing binaries and catalog false alarms as separately labeled
hard-negative cohorts in a later benchmark release. Do not count an eclipsing
binary transit as an empty/null light curve: a transit search should usually
detect its eclipse, while astrophysical classification/vetting is a separate
task. Future additions should include TOI dispositions with archived status
dates, pixel-level centroid labels, sector-held-out validation, and blind
injections over transit phase, depth, duration, gaps, and stellar variability.

# Data and provenance manifest

## Known-answer release

- Release directory: `data/benchmarks/generated-20261009/` (ignored by Git).
- Source catalogue: NASA Exoplanet Archive `pscomppars`, TESS-discovered
  confirmed transiting planets, TIC host identifiers; frozen CSV and query are
  in the release directory.
- Positive sample: 400 distinct host stars, one planet per host; 200
  calibration and 200 blind.
- Control sample: one median-reflected/inverted light curve per host; also 200
  calibration and 200 blind.
- Light curves: 3,824 MAST QLP sector products spanning Sectors 1–104, then
  stitched to one compact input curve per host. The metadata records product
  identifiers and SHA-256 hashes.
- Lightkurve: 2.6.0. The NASA catalogue snapshot hash and the builder/lockfile
  hashes are recorded in `provenance.json`.
- Inputs: 400 positive and 400 inverted curves; 800 rows in
  `truth-manifest.csv`. The label-free `run-manifest.csv` has 600 rows: the 200
  calibration null controls and all 400 blind cases. Calibration positives
  are excluded from this 600-case protocol and threshold selection. The audit
  also found 183 unreferenced outputs from earlier runs; those are excluded
  from all reported denominators.
- Search: 1 pass; `--threads 1`; 24 process workers; event grouping; period
  bounds 0.15–30 days; no explicit period count; default 250,000 period-trial
  cap.
- Output sidecars have one binary hash and one search configuration. All 600
  required input SHA-256 values match both truth metadata and sidecars; all 600
  candidate output SHA-256 values match sidecars; all 600 required candidate
  outputs exist.
- Candidate CSVs use the pre-detection-tier schema. Jackknife scores were
  reconstructed from the saved event-level S/N evidence using the same formula
  as the later source implementation. The binary SHA is recorded, but the
  executable itself and its source commit are not archived in this release.
- The audit found 183 extra, unreferenced result sidecars from earlier runs.
  They were not used: the runner and evaluator take their cases from the exact
  600-row run manifest.
- Total generated release footprint: about 2.4 GB. Keep it outside Git or
  upload to an immutable research-data archive before removing ignored files.

## Integrity report

The reusable audit command is:

```sh
.venv/bin/python scripts/audit_known_answer_benchmark.py \
  --run-manifest data/benchmarks/generated-20261009/run-manifest.csv \
  --truth-manifest data/benchmarks/generated-20261009/truth-manifest.csv \
  --provenance data/benchmarks/generated-20261009/provenance.json \
  --output data/benchmarks/generated-20261009/audit.json
```

The full machine-readable audit is
`data/benchmarks/generated-20261009/audit.json`. It verified exact run/truth
case coverage, no truth columns in the run manifest, uniform run settings,
matching data/output hashes, and `search_config_matches_provenance=true`.

## Limitations and source-data changes

The positives are confirmed-planet systems selected by prior discovery and
catalog inclusion. Their recovery is conditional on this sample, not survey
completeness. The inverted curves preserve target cadence and noise structure
approximately, but do not represent the population of eclipsing binaries,
background eclipsing binaries, blends, stellar variability, or other
astrophysical confounders.

The official MAST QLP documentation notes that from Sector 94 QLP uses the
TESS-Gaia Light Curve methodology. This release includes 276 product-sectors
from Sectors 94–104, so its input extraction is not fully homogeneous across
all years. Stratify or record extraction family/version when using this sample
for algorithm comparisons. See the [MAST QLP product page](https://archive.stsci.edu/hlsp/qlp).

No target-pixel or full-frame pixel data are included in the known-answer
benchmark. It evaluates light-curve detection and ephemeris recovery, not
centroid validation or source localization.

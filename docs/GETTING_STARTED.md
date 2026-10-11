# Astrosearch: getting started and running a survey

Astrosearch searches TESS light curves for **transit-like dips, periodic variability, and isolated events**. It is a search and triage system, not an object-confirmation system. Its primary input is a calibrated light curve; the Rust search does not localize a signal on the sky from pixels by itself. Every output is a candidate signal until follow-up rules out competing explanations.

## How a search works

1. The survey selects target TICs and TESS products from MAST (QLP and SPOC are separate reductions of TESS observations, not independent telescopes). It downloads sectors, preserves product/source metadata, quality, exposure, and sector information, and checks that usable measurements actually exist. Catalog coverage alone does not guarantee a usable light curve.
2. Each sector is normalized and checked, then usable sectors are searched with gaps kept as separate segments. The detector detrends and searches transit-like shapes, estimates correlated noise by segment, and performs five passes by masking the strongest event windows before searching for additional event patterns. Profile J reduces the period grid with cadence binning, a geometry-informed duration floor, a two-stage grid, and coarse fast-folding ranking.
3. Search scores rank signals. A tier label is meaningful only when the input product, binary, settings, and calibration contract match. A high score is not a probability that a signal is real. Candidate CSVs retain aliases and event-sharing alternatives for inspection.
4. Vetting checks known catalogs/TCEs, odd/even and secondary-event behavior, consistency across sectors, nearby/blended sources, pixels/centroids, and false-positive scenarios. Where possible, compare another reduction or algorithm. QLP versus SPOC is a useful reduction check, but both use TESS data. Pixel checks and TRICERATOPS/FPP require the relevant data and assumptions; missing evidence must be reported as missing.
5. A credible new-object claim needs reproducible measurements, source localization, an end-to-end calibrated search, and predictive independent follow-up. Publish or record predicted future event times before follow-up. Do not call a candidate confirmed from a search score alone.

The first-principles objective is **validated, reproducible signals per compute and follow-up budget**. Keep measured data provenance, detection sensitivity, source localization, and independent prediction separate. Do not tune thresholds on the blind evaluation split.

## Install and make a small first run

From the repository root:

```sh
uv sync
cargo build --release
```

The simplest workflow check uses the bundled example:

```sh
uv run python scripts/tess_survey.py \
  --targets data/survey/targets-example.csv \
  --output data/survey/first-run \
  --max-targets 1
```

This example is TIC 4206066, an existing unconfirmed benchmark signal. It checks installation, MAST access, product handling, and output shape; it is **not a fresh target, a discovery, or a truth standard**. The survey defaults to Profile J, QLP and SPOC, recent sectors (90+), a multi-year baseline, and five transit-search passes. Actual eligibility also depends on current MAST products and valid target metadata. Run a small pilot first and read its errors before increasing the target count.

To search a target list, provide a CSV with at least `tic_id,teff_k,distance_pc,radius_solar,tmag`; `stellar_density_g_cm3` is optional. Start with a bounded pilot:

```sh
uv run python scripts/tess_survey.py \
  --targets my_targets.csv \
  --output data/survey/my-first-survey \
  --max-targets 3
```

The default target filters favor nearby K/M dwarfs with long coverage. If preparing a list from target-yield rankings, join the ranked TICs back to the original TIC metadata; rankings alone may not include the stellar fields the survey requires. Check `--help` for adjustable sector, baseline, magnitude, temperature, distance, cadence, worker, author, and network settings.

The standalone Rust CLI also defaults to Profile J:

```sh
target/release/astrosearch --input lightcurve.csv --output candidates.csv
```

`--search-profile standard` selects the full-grid baseline. A fixed `--periods` grid requires `standard`. Profile J’s built-in detection threshold is calibrated only for its matching QLP five-pass contract; a result from another product or changed binary/settings needs matching calibration before the tier can be interpreted as calibrated. The survey records contract-match status. Do not treat an unmatched tier as validated evidence.

## Read results and decide what to inspect

For each run, inspect:

- `batch_summary.json`: per-target and per-author status, including unavailable, partial, and failed inputs.
- `run_metadata.json`: selected settings, calibration/contract information, counts, timing, and retry queues.
- `local_decisions.csv`: target selection/defer decisions and predicted yield; it is a prioritization heuristic, not a detection result.
- `TIC.../<author>/`: normalized sector curves, product manifest/checksums, stitched curve, candidate CSV, and search log.
- the shared SQLite survey ledger: which targets were completed, deferred, or are retryable. Keep it with the survey outputs for resumability and audit.

Only use data that reached the search stage. “No usable data” means product lookup/download/normalization yielded no usable measurements; “insufficient baseline” means the usable time span failed the configured requirement. Read the author-specific reason and logs; do not count these as successful searches or as null detections. MAST ping is not a useful health check: use the actual HTTPS/API requests and inspect their status and response.

For a candidate signal, first check the TOI/CTOI, SPOC TCE, and NASA Exoplanet Archive records. Then inspect its light curve, aliases, event count, odd/even depths, phase-0.5 secondary events, and per-sector depth/epoch consistency. Check Gaia neighbors, TIC contamination, and pixel centroids/localization; use pixel products where available. Run false-positive assessment with its inputs and assumptions recorded. Compare QLP/SPOC or TLS as an algorithm/reduction robustness check, not independent telescope confirmation. Record absent checks as absent. `candidate_preflight.py` and `ctoi_gate.py` are analysis/gating tools; they do not authorize or perform an ExoFOP submission.

Before trusting score tiers or yield rankings, use known-answer and injection/recovery runs with the exact production profile, full number of passes, data products, and blind split. Include inverted controls and realistic eclipsing-binary, blend, and systematic negatives. Report recovery and false-alarm rates with confidence intervals and paired gains/losses against the standard profile. An older calibration file is not automatically valid for a new binary or changed search contract.

## Run over time without repeating work

The survey shares a SQLite ledger and MAST cache by default, even when output folders differ. Each finished target is recorded so an interrupted run can resume; successful targets are normally deferred for 30 days, while failed targets can retry. The catalog row and search configuration are part of the ledger identity. Keep the same target catalog, output, ledger, and cache to continue the same campaign:

```sh
uv run python scripts/tess_survey.py \
  --targets my_targets.csv \
  --output data/survey/my-long-run
```

Use `--max-targets N` for bounded sessions, then rerun the same command later. Use `--revisit` only when you intend to reprocess completed targets now; `--revisit-after-days 0` disables the success cooldown. Use `--ledger` and `--cache-dir` to isolate campaigns. Avoid running two sessions against the same ledger unless deliberate: target leases prevent duplicate work, but a second process can report a target as busy. Worker defaults are conservative; raising target/download workers can saturate MAST, memory, or CPU without making searches faster. Increase them only after measuring end-to-end wall time and resource use.

For transient author/product failures, inspect the generated `author-retry-queue-*.csv` and rerun a bounded queue with the original settings. Keep the original run folder and ledger. Do not manually mark failed or partial work as complete. A restart should report what was retried, skipped, unavailable, timed out, or failed.

## Practical lessons and things to avoid

- Verify that actual time/flux/error/quality rows were downloaded. Metadata promises are not measurements; preserve sector gaps and product provenance.
- Start with a small pilot and inspect logs, usable baseline, thresholds, contract match, and runtime before a large survey.
- Profile J is optimized for throughput, but its duration floor encodes geometry assumptions and can miss very grazing or eccentric short events. Use standard or a purpose-built sensitivity analysis when those cases matter.
- Five-pass masking helps uncover additional event patterns but can mask overlapping signals. Compare alternatives and use injection recovery to measure the effect.
- Detrend each segment; do not let sector offsets masquerade as long-period variability. Re-measure transit depths with transits masked and fit a physical, finite-exposure model before interpreting radii.
- A detected dip can be an eclipsing binary, blend, instrumental artifact, or alias. Do not equate score, recurrence, or QLP/SPOC agreement with confirmation.
- Do not choose thresholds from the blind set, claim significance from a small convenience sample, or compare profiles at different false-alarm rates.
- Do not report counts, candidate parameters, runtime, or confidence intervals unless they were measured from saved outputs. Keep errors and timeouts in the report.
- Never call an unconfirmed search result a planet, never submit to ExoFOP, and never commit `.env` or credentials.

## Copy-ready prompts for an LLM operator

**First bounded search**

> Read `AGENTS.md`, `docs/GETTING_STARTED.md`, and the current CLI help. Check the repository state and existing survey ledger before running anything. Run a one-target Profile J pilot using `data/survey/targets-example.csv`, saving to a new output folder. Do not overwrite prior results. Report measured target/product counts, usable sectors and baseline, search status, calibration-contract match, runtime, errors/timeouts, and exact output paths. Treat all signals as candidates; do not call anything a planet or submit anything to ExoFOP. If prerequisites fail, report the failure and stop before scaling up.

**Resume a long survey**

> Continue the survey in `<output-folder>` using `<targets.csv>` and its existing SQLite ledger/cache. Read `run_metadata.json`, `batch_summary.json`, and the ledger first; do not repeat completed targets unless I explicitly request `--revisit`. Process at most `<N>` newly eligible targets this session. Use the same Profile J settings and calibration contract. Preserve author retry queues and report counts for searched, deferred, unavailable, partial, failed, retried, and timed-out targets from the saved run outputs. Do not invent measurements, call candidates planets, or submit to ExoFOP.

**Review one candidate**

> Review candidate `<id>` from `<run-folder>`. Follow the verification order in `docs/VERIFICATION_PROTOCOL.md`: provenance and usable data, known-catalog/TCE checks, aliases and multi-event support, odd/even and secondary checks, sector consistency, neighbors/blends and pixel localization, false-positive assessment, and prospective ephemeris prediction. State which checks were completed, their evidence/paths, and which are unavailable. Treat QLP/SPOC as reductions of the same telescope data. Separate observations from interpretation. Do not call it a planet, claim confirmation, or submit to ExoFOP.

**Evaluate a proposed algorithm change**

> Evaluate `<change>` against the current Profile J default from first principles. Preserve target/group-level train/calibration/blind separation. Use the full production pass count, matched source products, realistic positive injections added before preprocessing, and inverted plus eclipsing-binary/blend/systematic negatives. Compare standard and proposed methods at matched false-alarm rates with confidence intervals, paired gains and losses, runtime/CPU, and end-to-end recovery. Do not tune on the blind split. Save new calibration outputs separately; do not replace current calibration until the evidence and contract are reviewed. Report measured results only.

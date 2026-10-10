# astrosearch

`astrosearch` is a Rust candidate-triage tool for TESS-like stellar light curves. It reads one CSV/TSV light curve, performs shared robust trend estimation, and screens for periodic transit-like dips, sinusoidal variability, and isolated brightening/dimming events. All detectors reuse the ingested and preprocessed series, while each performs its own search over the data.

## Build and run

```sh
cargo build --release
./target/release/astrosearch --input lightcurve.csv --output candidates.csv
```

## TESS batch survey

The Python front end downloads public TESS light curves through Lightkurve/MAST, normalizes each sector independently, keeps pipeline reductions separate, and invokes the Rust search. Install Python dependencies with `uv sync`, then provide a TIC/CTL target-list CSV. Required column: `tic_id`. Recommended TIC fields are `teff_k`, `distance_pc`, `radius_solar`, and `tmag`; an optional `stellar_density_g_cm3` column enables a target-specific adaptive transit period grid. Targets without a usable temperature or distance are skipped. This workflow prioritizes cool, nearby, small stars with Sector 90+ observations and at least 1,095 days of total coverage. The recent sector can be combined with older sectors for a longer baseline. The priority score only sorts targets; it is not a detection probability.

```sh
uv sync
cargo build --release
python scripts/tess_survey.py --targets data/survey/targets.csv --output data/survey/run1
```

The example input is [targets-example.csv](/Users/danny/dev/astronomy/data/survey/targets-example.csv); it uses TIC 4206066, an unconfirmed candidate benchmark in `data/tic4206066`, so it validates workflow behavior and is not a fresh target or truth standard. Default authors are `QLP` and `SPOC`; results live in one folder per TIC and author, with one normalized CSV per sector, a manifest of product metadata/checksums, a stitched CSV, candidate CSV, and Rust search log. Light curves now retain `exposure_days` and `sector` columns so the refinement and injection tools can respect mixed cadence. These integrations use Lightkurve's MAST search and download interfaces; for large all-sky catalog queries, use MAST/astroquery to prepare a bounded TIC/CTL list first, then run the batch on a small pilot sample. QLP products use a 200-second minimum exposure by default; SPOC uses a separate 120-second minimum so standard 2-minute SPOC products are not silently excluded. Raise `--min-exptime-spoc` to reduce SPOC data volume. Each sector receives a robust adjacent-difference noise estimate normalized for exposure time; sectors above 5× that target's median are excluded from the combined search and recorded in the manifest. This is an outlier screen, not a universal threshold; inspect the ratios and tune `--max-sector-noise-ratio` (0 disables it). Tune `--min-sector`, `--min-baseline`, `--min-exptime`, and stellar cuts for the science question.

Before ranking a bounded TIC target list, attach real MAST coverage and restrict to targets meeting the survey's recent-sector and multiyear baseline requirements:

```sh
uv run astrosearch-enrich-targets targets-ranked.csv \
  --output targets-coverage.csv \
  --eligible-output targets-survey-eligible.csv \
  --cache data/survey/coverage-cache
uv run astrosearch-target-yield targets-survey-eligible.csv \
  --output targets-yield.csv --sky-cells sky-cells.csv
```

The coverage command caches per-target MAST results for resumability and records sector counts, observed days, and exact SPOC 2-minute sector IDs. Cached values are snapshots; delete a target's cache JSON to refresh it after MAST coverage changes. The survey still checks actual products and baseline before searching.

After a batch search, check all detection-tier rows against SPOC TCE files in one pass. The command downloads each shared MAST sector file once and records its URL and SHA-256:

```sh
uv run astrosearch-tce-preflight \
  --run-dir data/survey/run1 \
  --coverage targets-coverage.csv \
  --output data/survey/run1/spoc_tce_catalog_checks.json
```

By default this checks recent sectors (90+) with SPOC 2-minute coverage; no match is limited to the reported files and sectors and does not establish novelty.


Each successful Rust search is also added to `<output>/candidate_registry.sqlite` with the light-curve/candidate hashes, source/build identity, product manifest, and candidate measurements. Promotions remain human-evidenced: use `python scripts/candidate_registry.py promote ...` only after recording the relevant catalog, pixel, independent-pipeline, and follow-up evidence.

```sh
python scripts/tess_survey.py --targets targets.csv --output data/survey/pilot \
  --min-sector 90 --min-baseline 1095 --authors QLP SPOC --threads 24
```

The Rust box search is the fast screening stage. Refine its periodic-dip shortlist with a finite-exposure quadratic limb-darkened, circular-orbit model before spending time on manual vetting:

```sh
python scripts/refine_candidates.py \
  --lightcurve data/survey/pilot/TIC123456789/QLP/combined.csv \
  --candidates data/survey/pilot/TIC123456789/QLP/candidates.csv \
  --output data/survey/pilot/TIC123456789/QLP/refined.csv \
  --limb-darkening 0.4,0.2
```

The refiner reports model fit improvement over a per-event local baseline, fitted ephemeris and geometry, sector-by-sector transit depths, and a secondary-eclipse-at-phase-0.5 diagnostic. Add --mcmc-steps 500 after uv sync --extra inference for an optional emcee posterior sample. The default limb-darkening values are generic. Fit results are screening diagnostics only: the circular model, fixed limb darkening, and local least-squares covariance do not cover eccentricity, stellar-parameter uncertainty, correlated noise, aliases, dilution, or false-positive probability. Inspect MCMC chain mixing and use stellar priors before reporting posterior uncertainties.

Each full Rust search can be substantial on a multi-year light curve. Start with a handful of targets, review the period-grid warnings in `search.log`, and expand only after measuring runtime and recovery for that target class.

For an independent limb-darkened Transit Least Squares search, install the optional dependency and use the iterative screen, which masks each strongest transit before searching again:

```sh
uv sync --extra tls
python scripts/tls_search.py --lightcurve data/survey/run1/TIC123456789/QLP/combined.csv \
  --output data/survey/run1/TIC123456789/QLP/tls-candidates.csv \
  --min-period 0.15 --max-period 30 --max-planets 5
```

Use tess_survey.py --tls to run TLS per target. TLS false-alarm estimates assume white noise; compare with --invert controls and empirical nulls. This independent detection path does not validate candidates.

For a frozen known-planet comparison, run TLS on the same calibration controls
and blind light curves, calibrate its SDE threshold on calibration controls,
and optionally compare its paired blind outcomes with a saved astrosearch
evaluation report:

```sh
python scripts/run_tls_benchmark.py \
  --truth-manifest data/benchmarks/generated-20261009/truth-manifest.csv \
  --run-manifest data/benchmarks/generated-20261009/run-manifest.csv \
  --output data/benchmarks/tls-20261010 \
  --reference-report data/benchmarks/generated-20261009/evaluation-jackknife.json
```

This runs one TLS analysis per required case and can take hours. The old blind
sample has already been inspected, so the result is a retrospective method
comparison only and must not be used to tune production settings.

Sector 90+ is a recency filter, not evidence that targets or signals have not been searched. QLP processes FFI light curves and performs multi-sector searches; prioritize documented coverage gaps or measured residual sensitivity, and audit current catalogs before calling a signal novel.

### Resuming surveys across sessions

The batch runner keeps a shared SQLite ledger at `data/survey/survey-ledger.sqlite` and a shared MAST download cache at `data/survey/mast-cache`, even when each run writes to a different `--output` folder. Target results are committed after each target finishes, so an interrupted run can resume without treating a half-finished target as complete. By default, successfully processed targets are deferred for 30 days; failed targets can be retried immediately. A changed catalog row or search configuration makes the target eligible again. After the cooldown, the runner checks the current MAST products and can pick up newer sectors. Use `--revisit` to refresh targets immediately, `--revisit-after-days 0` to disable the cooldown, or pass `--ledger` and `--cache-dir` to keep independent surveys separate. For example:

```sh
python scripts/tess_survey.py --targets targets.csv --output data/survey/session-2 \
  --max-targets 25
python scripts/tess_survey.py --targets targets.csv --output data/survey/refresh \
  --revisit
```

The ledger is a work tracker, not a claim that a target or signal is novel. Keep the ledger file with the survey outputs so it can be resumed or audited later.

Only one batch process sharing the same ledger can hold a live claim for a TIC at a time. An abandoned claim is recoverable after 24 hours; Ctrl-C and SIGTERM mark in-progress targets retryable immediately. MAST HTTP requests use a 300-second connect/read timeout by default and bounded retries for connection errors, HTTP 429, and server errors; use `--network-timeout SECONDS` to change the request timeout. Rust searches have a 2-hour per-target timeout by default; use `--search-timeout SECONDS` to change it. Target failures produce a partial/failed run status and a nonzero exit code. Inspect `batch_summary.json`, `run_metadata.json`, and the ledger before treating a batch as complete. Batch CSVs and manifests are written to temporary files and atomically published; manifests include hashes for source sectors and search outputs to support audit and resume checks.

By default, eligible targets are ordered by the **expected number of new planets this search can detect** around each star (`--ranking yield`, implemented in `scripts/target_yield.py`, also available as `astrosearch-target-yield targets.csv --output ranked.csv [--sky-cells cells.csv]`). For each star the model integrates over orbital period (the survey's `--min-period` to `--max-period`) and planet radius (0.5–20 R⊕):

- **Planet occurrence:** the SAG13 occurrence rates (Kopparapu et al. 2018), reduced inside about 10 days (Mulders et al. 2015). Small planets are twice as common around M dwarfs as around Sun-like stars; giant planets are half as common.
- **Transit probability:** the chance that the orbit is aligned to transit. For stars that already have a known transiting planet (`known_planet_periods`), this becomes the much higher probability for a nearly coplanar system.
- **Detection completeness:** the measured injection-recovery completeness of the deployed detection tier (`data/calibration/completeness.json`).
- **Discount for earlier searches:** planets that SPOC or QLP probably already found are discounted (`data/calibration/discovery_priors.json`). This step uses cited assumptions, not measurements.

Predicted signal-to-noise uses the 1-hour noise measured against TESS magnitude on 384 QLP light curves (`data/calibration/noise_vs_tmag.json`, built by `scripts/calibrate_noise_model.py`), or `cdpp_1h_ppm` when given. Dilution is reduced by `1 + contamination_ratio`. Observing time comes from `observed_days` or `sectors` when present, otherwise from `--assumed-sectors` (4).

The result is an expected count from population statistics. It says nothing about whether any particular star has a planet. Typical values are 10⁻⁴ to 10⁻² per star.

`--ranking snr` restores the older ranking, a relative fixed-radius-planet signal-to-noise proxy (CDPP, or a TESS-magnitude photon-noise fallback).

`--search-profile fast` runs grid-study variant J, which combines four settings:

- 10-minute inverse-variance binning;
- physical transit-duration limits;
- a 7× coarse period grid, ranked by fast folding and refined around the best 512 cells;
- a cap of 2 million period trials.

On 158 known-answer cases (78 planets, 80 noise-only), J used 7× less CPU than `standard`: 2,776 s against 20,563 s, and 6 minutes against 110 minutes of wall time. At matched false-alarm rates it recovered every planet that `standard` found, plus 3 more at 5% and 4 more at 10%. That gain is not statistically significant (p = 0.25 and 0.13). `standard` stays the default until injection completeness has been rechecked under `fast`. `--periods` cannot be combined with `fast`.

`local_decisions.csv` records the yield for each eligible target, along with whether the target was selected, deferred by the budget, or suppressed by the ledger. Both rankings are deterministic and need no external service.

## CTOI readiness gate

`scripts/ctoi_gate.py` runs one transit candidate through the checks an exoplanet vetting team would expect before a Community TESS Object of Interest (CTOI) submission, and writes `gate.json`, a readable `gate.md`, a draft `ctoi_row.csv` with the ExoFOP upload fields, and `ephemeris_windows.csv`. It never submits anything and never calls a candidate a planet.

```sh
uv sync --extra inference --extra fpp --extra pointing   # batman, emcee, triceratops, tess-point
cargo build --release
./target/release/astrosearch --input data/tic4206066/combined.csv --output cands.csv
python scripts/ctoi_gate.py --tic 4206066 --lightcurve data/tic4206066/combined.csv \
  --candidates cands.csv --rank 21 --products qlp_s32.csv \
  --pixel-reports pixels/s32/pixel_vet.json --out-dir gate/tic4206066
```

| # | Criterion | Automated checks |
| --- | --- | --- |
| 1 | Signal in 2+ sectors and 2+ independent products | events in at least two sectors; per-sector depths consistent (χ² p ≥ 0.01); two pipelines (from the light-curve parts or `--products`) agree |
| 2 | Physical fit, consistent star | limb-darkened transit fit preferred by ΔBIC ≥ 10; transit-implied density 3π(a/R*)³/(GP²) within 3σ of TIC (MCMC posterior uncertainties by default) |
| 3 | No false-positive indicators | secondary eclipse searched at every phase (model-shift style: SNR scaled by the scan's own red-noise factor, a look-elsewhere threshold for ~P/T phases, a single-event jackknife, and a margin over the strongest positive excursion; 3σ at phase 0.5); odd/even depths within 3σ; not grazing (b < 1 − Rp/R*); events not dominated by momentum dumps or scattered light; no difference-image centroid offset that is both beyond one TESS pixel and ≥ 3σ (run automatically on up to 2 SPOC sectors, target position propagated with proper motion); no TIC neighbor within 2 pixels (separations at the TESS epoch) both bright enough to produce the depth and within the Vanderburg et al. 2019 ingress/egress magnitude limit; Gaia DR3 RUWE < 1.4 (matched by source_id) |
| 4 | Not already known | no ephemeris match in the TOI, CTOI (ExoFOP), confirmed-planet, or requested SPOC TCE tables |
| 5 | Significant and not a false positive | jackknife score at or above the calibrated detection tier (`data/calibration/detection_threshold.json`); empirical per-target FAP 95% upper bound < 1% from inverted and block-shuffled null light curves (300 by default); TRICERATOPS FPP < 1% |
| 6 | Testable prediction | ephemeris with uncertainties and predicted windows; publishing them with a timestamp and submitting the CTOI are human actions |

Before the diagnostics, other transit candidates scoring at least half the candidate's score are masked unless their period is commensurate (1/3, 1/2, 1, 2, 3×) with it, so a sibling planet cannot pose as a secondary eclipse or odd/even difference while a same-period secondary stays visible (`--sibling-score-fraction`, 0 disables). Full BJD inputs are normalized to BTJD. `gate.json` records the git commit, input SHA-256 digests and package versions, and `gate.md` ends with the [verification protocol](docs/VERIFICATION_PROTOCOL.md) checklist.

Each criterion is PASS, FAIL, NOT_EVALUATED (missing data, tool, or network), or ACTION_REQUIRED. The verdict is "ready for human review" only when every automated criterion passes. Depths from products without a pipeline crowding correction are scaled by 1 + the TIC contamination ratio before deriving the planet radius, insolation, and equilibrium temperature (Monte Carlo over TIC uncertainties; circular orbit; albedo 0.3). Two-transit signals are reported by `astrosearch` as `duo_transit`; the gate enumerates their allowed period aliases (`duo_aliases`), ruling out aliases that predict a transit in data showing no dip. Statistical validation additionally needs FPP < 0.015, NFPP < 0.001, and high-resolution imaging, which this tool does not provide.

Supporting tools:

- `scripts/export_tess.py` now writes a `<file>.meta.json` sidecar per sector (product, flux column, whether the pipeline corrected crowding, momentum-dump and scattered-light times); the gate reads these.
- `astrosearch --output X` writes `X.provenance.json` with the tool version, git commit and dirty flag, the full configuration and its SHA-256, and the input file's SHA-256.
- The candidate CSV's `jackknife_score` column is the score with the strongest single event removed, and `tier` labels transit candidates `detection` (jackknife score ≥ `--detection-threshold`, default 23.3) or `screening`. Score 6 (`--min-transit-score`) is only a screening floor: 97% of the known-answer benchmark's inverted null controls exceed it, and their best candidates are usually carried by one event. On the blind split, at thresholds fixed on the calibration split, the jackknife tier recovers 86/200 planets with 14/200 null triggers, against 61/200 and 8/200 for a raw-score tier; at equal blind false-alarm rates it recovers more at every rate tested (2–10%). `scripts/calibrate_detection_threshold.py` recomputes the tier from the benchmark.
- `scripts/evaluate_secondary_scan.py`, `scripts/injection_recovery.py`, `scripts/baseline_bls_benchmark.py` and `scripts/compare_fap_passes.py` reproduce the review-board experiments (see [docs/REVIEW_IMPROVEMENTS.md](docs/REVIEW_IMPROVEMENTS.md)).
- `scripts/gp_detrend.py` flattens active stars with a celerite2 Gaussian process (transits optionally masked); search its output with `--trend-method none`, and compare with the biweight default using injection-recovery before preferring it.
- `scripts/benchmark_known_systems.py build|run` checks recovery and vetting on known planets (LHS 3844 b, L 98-59 b–d, TOI-700 b–c) and deep TFOPWG false-positive TOIs, with ephemerides pulled from the NASA Exoplanet Archive into `data/benchmarks/known_systems.csv`.

## Candidate preflight

For a candidate file, generate a dated report that checks the NASA Exoplanet Archive TOI and confirmed-planet tables by TIC ID, compares catalog ephemerides with period/phase aliases, and creates a checklist for remaining checks:

```sh
python scripts/candidate_preflight.py --tic-id 123456789 \
  --candidates data/survey/run1/TIC123456789/QLP/candidates.csv \
  --tce-sectors 90,91,92 --output data/survey/run1/TIC123456789/preflight.md
```

The preflight fetches SPOC TCE CSVs for the supplied sectors and links to ExoFOP. A missing TCE row only means no TIC match in the files that were fetched; also inspect available multi-sector coverage. An empty TOI query does not check the CTOI list or establish novelty. Independently compare QLP and SPOC; fit odd/even depths, secondary eclipse, and sector-by-sector depths; inspect Gaia neighbors/RUWE, TIC contamination, and target-pixel centroids; then run TRICERATOPS with documented field-star, aperture, and stellar inputs. Do not interpret the astrosearch score as an FPP. Only publish predicted transit windows after fitting ephemeris uncertainties, and submit a CTOI after assembling the validation package.

For a per-sector pixel localization diagnostic, compare in-transit and out-of-transit target-pixel images:

```sh
python scripts/pixel_difference_check.py --tic-id 123456789 --sector 98 \
  --period 3.18 --epoch-btjd 1468.42 --duration-days 0.12 \
  --ra-deg 123.456 --dec-deg -12.345 \
  --output data/survey/run1/TIC123456789/pixel-vet-sector98
```

Inspect the difference and S/N images manually; a centroid offset is a warning, not a disposition. The tool does not replace SPOC/DV difference-image products, Gaia-neighbor analysis, or TRICERATOPS. Run it independently for each sector containing transits.

The reusable known-answer benchmark builder now freezes confirmed TESS planet light curves and paired inverted controls, then calibrates a score threshold on a separate host-star split. The repository includes the workflow but no downloaded benchmark release yet; do not treat example evaluator output as validated completeness or a calibrated false-alarm rate. See [data/benchmarks/README.md](/Users/danny/dev/astronomy/data/benchmarks/README.md) for generation, blinded execution, evaluation, and statistical limitations.

Once a transit model provides ephemeris errors, prepare a BJD_TDB observing schedule. Astrosearch's screening epoch alone has no fitted uncertainty:

```sh
python scripts/predict_transits.py --period 3.18 --epoch 1468.42 \
  --period-sigma 0.0001 --epoch-sigma 0.01 --after-btjd 5000 \
  --through-btjd 5100 --output predictions.csv
```

The example values above are placeholders; use the transit-model fit and include covariance when available. The generated windows are planning material, not a public pre-registration. Publishing and CTOI submission remain deliberate researcher actions.

The CLI also provides a reflected-flux systematics control: `astrosearch --input lightcurve.csv --output inverted-candidates.csv --invert`. Compare it with the original search; it is an empirical null check, not an astrophysical false-positive model. The real-curve validation workflow already runs this inversion and leave-one-season-out ephemeris prediction.

The injection/recovery example in `examples/mixed_cadence_injection.rs` reports completeness and null-curve alarm fractions over several score thresholds. It supports white noise (`red-rms-fraction 0`) and a toy within-orbit AR(1) component (default fraction 0.5, correlation 0.8). For example, `cargo run --release --example mixed_cadence_injection -- 24 12 20 0.5` runs 12 injections per depth and 20 null curves. This is a smoke-scale sensitivity study; repeat with larger samples, multiple correlation times, and real light-curve nulls before using it to calibrate a survey-wide threshold.

### Validation on real light curves

Run randomized injection/recovery, an inverted-flux null search, per-segment circular-shift null searches, and a leave-one-observing-season-out prediction check on an extracted target:

```sh
python scripts/real_curve_validation.py \
  --sector-dir data/survey/run1/TIC123456789/QLP \
  --output data/survey/run1/TIC123456789/QLP/validation \
  --min-period 1 --max-period 15 --injections-per-depth 50 --null-trials 100 --threads 24
```

The script writes per-sector noise diagnostics, the filtered combined curve, real and inverted candidate lists, time-shifted null maxima, injection recovery outcomes, holdout-season results, and a JSON summary. Injections randomize period, phase, duration, depth, and impact parameter. The default is a numerically integrated quadratic limb-darkened transit, with seven-point finite-exposure integration; trapezoid and box controls are available with `--injection-shape`. The default limb-darkening coefficients are generic placeholders (`u1=0.4,u2=0.2`). The new survey exports include actual exposure duration; older products fall back to median timestamp spacing. Supply stellar-model TESS-band coefficients for stronger physical fidelity. The geometry assumes a circular orbit with linear projected motion. Recovery is recorded at any candidate rank and summarized separately for top 1, top 5, and top 10; Wilson intervals make the uncertainty from small samples visible. Results also report recovery among injections with at least two epochs containing two or more in-transit samples, so data-gap-limited trials are not conflated with search failures. Circular-shift nulls independently rotate flux and uncertainties within each continuous observing segment, preserving local correlation approximately while breaking cross-segment phase coherence. Their plus-one tail fraction reports its finite-sample resolution and is an exploratory maximum-score comparison, **not** a calibrated false-alarm probability. Each holdout fold discovers candidates using all other observing seasons, then measures predicted transit windows only in the withheld season. These are target-specific pipeline diagnostics, not survey completeness or planet validation. Use at least 50 injections per parameter bin as an initial diagnostic, then increase the number across many targets and parameter bins before making population-level claims.

To compare median and Tukey-biweight trend estimates with paired injection phases and null trials:

```sh
python scripts/compare_detrenders.py \
  --sector-dir data/survey/run1/TIC123456789/QLP \
  --output data/survey/run1/TIC123456789/QLP/detrender-comparison \
  --injections-per-depth 50 --null-trials 100 --threads 24
```

The Rust default is a time-windowed Tukey biweight trend with a 1.05-day window, three times the maximum searched transit duration. --trend-days changes the filter width; --segment-gap-days independently controls gap boundaries. Every detected candidate is re-detrended with its transit windows excluded, and its depth is remeasured as the median of nearby out-of-transit samples (1.5–3 durations from each event) minus the median in-transit value. The remeasured depth replaces the box-fit depth only when at least 2 in-transit and 4 nearby samples exist and the result is positive; those rows say "depth remeasured after transit-masked detrending" in `note`, and other rows keep the box-fit depth that `score` uses. Keep the median method available for paired recovery comparisons because no detrending choice is optimal for all target variability and transit durations.

Requires Rust 1.87 or newer. Input must have a header with `time` (or `btjd`/`bjd`) and a positive `flux`. Optional columns: `error`/`flux_err` and `quality` (only zero quality is retained). Times are in days and errors must share the flux units. Invalid/nonpositive flux rows are skipped and counted; unavailable or invalid errors use a scatter estimate. Already mean-subtracted light curves must first be converted to positive relative flux. Quoted fields, escaped quotes, BOM headers, CSV, and TSV are supported; multiline quoted records produce an explicit error. If both SAP and PDCSAP columns exist, PDCSAP is preferred and its matching uncertainty column is used. Explicit `flux` takes precedence over both. The CLI rejects an output path that refers to the input file, including symlinks and Unix hard links.

```sh
astrosearch --input target.csv --min-period 0.2 --max-period 40 --periods 4000
```

Both periodic scans automatically use available CPU parallelism. Use `--threads 1` for serial execution, `--threads 16` to limit workers, or `--threads 0` (the default) for automatic selection. Each worker keeps private phase buffers while sharing the light curve; dynamic chunks balance varying period costs and fast/slow CPU cores. Trial results are restored to serial order before ranking, so scheduling does not change candidates. The two periodic detectors run sequentially, avoiding two full worker pools competing at once. If an external batch runner processes several stars concurrently, limit threads per process to avoid oversubscription. Automatic availability is an OS estimate and may vary under resource limits ([Rust API](https://doc.rust-lang.org/std/thread/fn.available_parallelism.html)).

## Search design

Exact duplicate timestamps are median-combined before the minimum-observation check. High-amplitude samples are preserved because they may be real flares or dips.

The light curve is first split into independent segments, such as TESS sectors or orbits observed with one data product. A segment ends at a gap longer than the larger of segment-gap-days and five of its own cadences, or where the sampling cadence changes by more than 1.5x and the new spacing persists for 50 samples; short dropouts and isolated samples do not start a segment. Each segment has its own median cadence. The CLI prints each segment time range, observation count, cadence, and noise to stderr, and Report::segments returns the same summary.

Trend anchors are medians of occupied bins at their observed time centers, computed separately within each segment; interpolation is linear and never crosses a segment boundary. No empty bins are allocated across gaps. A coarser trend feeds the variability search. Flux errors are divided by the same fitted trend as the flux. Within each segment, neighbor differences (excluding gaps over three of that segment's cadences) give a robust local-scatter estimate and a conservative floor on quoted uncertainties. A segment with fewer than 20 differences borrows the pooled estimate of segments within 1.5× of its cadence, then the whole light curve. Per-segment noise matters for combined multi-year TESS data: a single global floor estimated from 200-second samples cut the weight of a 30-minute sector of TIC 4206066 by 11×. A separate duration-binned beta estimate measures correlated noise in the detrended transit series; neither estimate captures all trend uncertainty or nonstationary systematics.

The variability frequency ceiling uses the global median cadence. The transit scan accounts for the duration floor described below, while mixed-cadence inputs still use that global cadence as their minimum sampling scale.

The transit search uses a weighted, phase-folded box scan with a clipped out-of-window background and an approximate two-level SNR; background values are clipped while in-window transit values are kept intact. After finding the strongest signal, the search masks its observed event windows and repeats up to --transit-search-passes (default 5), then groups aliases by shared events and labels integer/half-integer harmonics in alias_ephemerides. Each pass contributes its strongest candidate, and the remaining `--transit-candidates` slots come from the last pass's ranked list. Masked samples are excluded from later passes (zero weight, not counted as in-window points) rather than replaced with a flat baseline, so they neither dilute nor fake-support later windows, and the correlated-noise β profile is re-estimated from the unmasked samples after each pass. Every pass reruns the full period search, so runtime grows almost linearly with passes; correlated-noise histograms and per-event bookkeeping are built only for trial periods with a window above `--min-transit-score`, which roughly halves each pass. On TIC 4206066 (quiet machine, 24 cores), the 5-pass default takes 20 s for sector 32 (17,094 samples) and 65 s for the combined three-sector curve. This is a sequential multi-planet screen; masking can still hide a genuinely overlapping transit, so compare with unmasked and TLS results. The period grid limits accumulated timing drift to half the shortest searched duration when uncapped. Without target metadata it retains the prior uniform log-period grid. With --stellar-density RHO (mean stellar density in g/cm^3), it samples uniformly in a duty-cycle coordinate whose local period spacing grows with the minimum geometric transit duration. The duration floor assumes impact parameter <= 0.95, eccentricity <= 0.9, and a radiusless planet approximation; these are explicit search priors, not universal limits. Omit the density option when grazing or more eccentric transits must remain in scope. The target-specific floor is never below two median cadences. Transit periods are limited to at least three times the duration floor (six cadences in the default mode) and at most the observation span. Actual observed epochs decide whether repetition exists, so two events can be found even when their period exceeds half the baseline. The largest searched duration is the larger of the duration floor and the smaller of 12% of the period and 0.35 days. The CLI reports total transit trials across all masked passes.

Weak events can accumulate evidence: an accepted window needs aggregate score at least 6, at least two dimming epochs with two samples each, and dimming in at least 60% of covered epochs. No individual 3-sigma requirement is imposed on every epoch. To screen out isolated events, the strongest epoch must carry less than 80% of the positive signal contribution; the effective number of contributing epochs must be at least the larger of 1.5 and half the covered epochs; and removing the strongest epoch must leave a combined score of at least 3. These are uncalibrated screening heuristics. Sparse coverage and unequal epoch sampling can reduce completeness.

Phase buffers and weighted residuals are reused across periods; circular prefix sums score windows without revisiting each point. The default grouping mode groups ephemerides that explain at least 80% of the same observed transit events. The highest-scoring representative remains a candidate; distinct alternatives are preserved in alias_ephemerides with harmonic labels rather than consuming candidate slots. This is a grouping heuristic, not proof that two solutions come from one astrophysical source. Up to 16 alternatives are retained per group. The sinusoid scan uses a frequency recurrence in blocks instead of evaluating sine and cosine for every sample at every trial, and subtracts each segment weighted mean before scoring to suppress sector-offset artifacts.

Periodic-dip CSV rows include `event_evidence` entries encoded as `midpoint:depth_ppm:independent_snr:point_count`, separated by semicolons. They describe observed support for that candidate; the per-event S/N assumes independent errors and is diagnostic only. `score` and `white_score` retain their existing meanings in both grouping modes. Candidate ordering remains by the selected search score; event grouping changes candidate slots and keeps alternative ephemerides attached to the representative.

Compare the legacy and event grouping modes on identical seeded calibration cases:

```sh
cargo run --release --example correlated_noise_calibration -- --scale 0.02 --real-limit 4 --grouping ephemeris --out target/calibration-ephemeris
cargo run --release --example correlated_noise_calibration -- --scale 0.02 --real-limit 4 --grouping events --out target/calibration-events
```

The harness uses the same fixed seeds in both runs and includes white/correlated noise, glitches, two-event planets, two-planet systems, and inverted real curves when available. Use `--real-limit N` to bound the number of sorted real-curve controls in a pilot. Reports include top-1/top-5/top-12 recovery at the 5% null-derived threshold, score-threshold recovery, Wilson intervals, runtime, and raw per-case outcomes; aliases count as matching ephemerides but use their parent candidate-group rank. The scaled run is exploratory; use the full sample and inspect per-star null behavior before drawing survey-level conclusions. Keep the legacy mode available until paired recovery, false-alarm, alias, and multi-planet results show that event grouping is safe for the target population.

The weighted floating-mean sinusoid search respects both period bounds, uses a median-cadence frequency ceiling, and samples at least five frequencies per inverse-baseline resolution element. Weights are cached, and nearly singular projections are rejected. A heteroscedastic event scan includes both same-sign 2.5-sigma shoulders around a 4.5-sigma peak, requires `min-event-points`, and splits at segment boundaries and at gaps over three cadences of the event's own segment, so consecutive 30-minute samples can form an event inside a light curve dominated by 2-minute data. Only isolated dips contained within a strong periodic window are suppressed as duplicates.

Periodic-dip candidates include an `odd_even_sigma` diagnostic comparing aggregate depths in alternating observed cycles. A value at or above 3 is flagged as a possible eclipsing binary. An unavailable diagnostic is blank, rather than reported as zero. Correlated scoring inflates the in-window variance using segment-specific beta estimates measured at transit-like timescales; white scoring remains available for comparison. This approximate correction does not model timing variations, all out-of-window/systematic covariance, or missing-cycle aliases, and scores still need target-specific null calibration. Both periodic searches warn when their required grid exceeds the 250,000-trial cap, since narrow peaks may then be missed. `periods` is a minimum requested grid size (10–250,000); `durations` accepts 1–128, and `min-event-points` accepts 2–1,024. Floating-point options must be finite and positive.

This first version searches **one light curve at a time**. It can screen planet-like transits, eclipsing-binary-like dips, stellar/asteroid rotation if an asteroid light curve has already been extracted, periodic stellar variability, flares, and isolated dimming. It does not read TESS FITS pixels, find moving brown dwarfs in WISE images, detect image-level transients, or distinguish exocomets from other asymmetric/isolated dips. Those require separate image or event-shape modules. The program is a proof-of-concept, not yet a validated survey pipeline.

## Interpretation

Scores are ranking statistics, not calibrated probabilities or confirmation. Different detector classes have different templates and trial counts, so `rank` is only output ordering and `class_rank` is rank within a detector. Neither is a catalog of unique discovered objects. `epoch` is the center of a transit window or isolated excursion in the input time system; a sinusoid's epoch is its analysis time origin. Event duration includes one cadence of the event's segment for sample coverage. These center conventions replace the earlier start-time convention.

The correlated score uses an approximate duration-dependent noise correction, but neither score corrects for trend uncertainty or the number of stars/periods searched. Median detrending can attenuate long signals. Transit-like signals can be eclipsing binaries, blends, instrumental artifacts, or stellar variability. Before claiming a candidate, compare pixel-level data, known-object/candidate catalogs, independent sectors and pipelines, and follow-up observations. The auxiliary real-curve validation script now provides finite-exposure, quadratic limb-darkened injections, rank-aware recovery counts with Wilson intervals, coverage-conditioned recovery, inverted-flux controls, and circular-shift nulls. It still does not provide a survey-wide completeness surface, a calibrated false-alarm probability, formal Bayesian model comparison, or pixel-level vetting. TESS products are FITS; convert a selected light curve to CSV for this prototype. Pixel-level centroid and difference-image checks are outside its input model, so the odd/even flag is not a substitute for SPOC/DV or pixel-level vetting.

## Library and checks

The CLI is a thin wrapper around `astrosearch::analyze_file` and `astrosearch::analyze`. `analyze` takes a vector of `Point { time, flux, error }` and a `Config`, returning a `Report` with candidates, observation count, baseline, per-segment cadence/noise summaries, and warnings. `Candidate::kind` is a `CandidateKind` and alias labels are `AliasType` values; both print the same strings as the CSV and compare equal to them (`kind == "periodic_dip"`).

| Module | Responsibility |
| --- | --- |
| `config`, `types` | Options with defaults and validation; samples, candidates, and reports |
| `input` | CSV/TSV parsing, quality filtering, duplicate merging |
| `segments`, `preprocess`, `rednoise` | Cadence/gap segments, robust trends, per-segment noise and weights, duration-dependent β |
| `grid`, `phase`, `epochs`, `transit` | Period grid and duration floor; phase-folding buffers; per-event statistics and vetting; the multi-pass box search |
| `grouping` | Ephemeris deduplication, event-overlap grouping, alias labels |
| `variability`, `single_events` | Sinusoid search; isolated dips and flares |
| `refine`, `output` | Depth remeasurement after transit-masked detrending; CSV rendering |
| `parallel`, `math` | Deterministic parallel trial scheduling; medians and robust scatter |

```sh
uv sync --frozen --extra inference          # locked environment; installs this project editable
cargo fmt --check
cargo clippy --all-targets -- -D warnings
cargo test --release                        # also builds target/release/astrosearch
uv run --frozen python -m coverage run -m unittest discover -s tests -p 'test_*.py'
uv run --frozen python -m coverage combine && uv run --frozen python -m coverage report --fail-under=60
```

GitHub Actions runs exactly these steps on pushes and pull requests.

The Python side is an installable package: `scripts/vetting` plus the operational tools as modules. `uv sync` installs it editable. It is editable because the gate runs `refine_candidates.py` and the release binary from the checkout. Console commands: `astrosearch-gate`, `astrosearch-survey`, `astrosearch-export`, `astrosearch-pixel-check`, `astrosearch-refine`, `astrosearch-preflight`, `astrosearch-registry`, `astrosearch-calibrate-tier`, `astrosearch-benchmark-run`, `astrosearch-benchmark-evaluate`, `astrosearch-known-systems`, `astrosearch-enrich-targets` and `astrosearch-tce-preflight`. One-off research scripts stay in `scripts/` outside the package and outside the coverage gate. The gate is 60% of package lines; 67% at the time of writing, with `tess_survey` at 85%. The Python suite includes tests that need the release binary (set `ASTROSEARCH_BIN` to use another build; they skip if it is missing):

- `test_parser_parity.py`: differential fuzzing of the Rust and Python light-curve readers on randomized messy files.
- `test_cli_fuzz.py`: CLI mutation fuzzing.
- `test_gate_integration.py`: an end-to-end search → offline gate run.

Rust property tests (`tests/properties.rs`) check invariance to input row order and to invalid rows, plus the CSV shape. Every networked script installs finite request timeouts (`scripts/vetting/network.py`); the gate's is set with `--network-timeout`.

Regression checks cover shallow repeated transits, two nearby planet periods, two transits separated by more than half the baseline, null/isolated-event controls, flux/time invariance, CSV parsing and flux/error pairing, sparse gaps, invalid configurations and duplicate-only data, error normalization, event shoulders, source-file overwrite protection, and mixed-cadence segmentation (sector gaps, cadence switches in both directions, dropouts, tiny-segment noise fallback, per-segment noise and weights, event runs across segment boundaries, shallow-transit recovery, noise-only nulls, and serial/parallel bit identity). Passing these checks demonstrates the specified behaviors; it does not establish survey completeness or reliability.

## Mixed-cadence validation

> Measured with the defaults of the time: 0.5-day median trend, one transit search pass, and ephemeris deduplication. To reproduce these numbers with the current build, add `--trend-method median --trend-days 0.5 --transit-search-passes 1 --candidate-grouping ephemeris`.

`data/tic4206066/` holds real TESS light curves of TIC 4206066 at three cadences: sector 6 at 1800 s, sector 32 at 120 s, and sector 98 at 200 s, plus their time-sorted combination. The data README gives provenance and checksums. A transit candidate at P = 3.18 d (~500 ppm, unconfirmed) has been reported for this star, so its rank and score are a known-answer check. `scripts/review_runs.sh` reruns everything below. "Before" is the build immediately preceding per-segment noise, which already had parallel scans. Timings are single runs on the Mac Studio below with automatic threads.

Use `python scripts/known_answer_benchmark.py --reference-epoch PUBLISHED_BTJD_T0` as a search-recovery release gate, replacing the placeholder with a published transit midpoint in the input BTJD time system. It requires a periodic dip whose accumulated period drift across the baseline and phase offset from that reference epoch are each within one hour, with at least two events. The benchmark deliberately requires an epoch as well as a period: period agreement alone can match an unrelated dip train. No reference epoch is stored in this repository, so obtain it from the source paper or its transit-time table rather than guessing. The current default build is expected to fail this gate, as documented below. Treat that failure as an unresolved sensitivity problem; do not weaken the gate or change the benchmark ephemeris merely to make it pass. The report is written beside the candidate CSV. This checks signal recovery, not planet validation.

```sh
cargo build --release
sh scripts/review_runs.sh                      # full, narrow, inverted-null, serial/parallel runs
cargo run --release --example mixed_cadence_injection -- 0 12
```

Segmentation of `combined.csv` finds eight segments, one per TESS orbit: two at 1800 s (202–204 ppm), two at 120 s (843–850 ppm), and four at 200 s (598–681 ppm). Before this change, one 734 ppm floor, estimated only from 200-second neighbors, applied to every sample. That cut sector 6's weight by about 11×.

| Run (P = 3.0–3.4 d unless noted) | Before score / depth | After score / depth |
| --- | ---: | ---: |
| Combined, `--trend-days 0.5` | 13.06 / 354 ppm | **15.46 / 367 ppm** |
| Combined, `--trend-days 1.0` | 13.62 / 737 ppm (23 epochs) | **23.32 / 554 ppm (28 epochs)** |
| Inverted combined (null), trend 0.5 | best dip 7.50 | best dip 7.77 |
| Sector 6 / 32 / 98 alone, trend 0.5 | 8.88 / 9.67 / 9.67 | 8.88 / 9.63 / 9.67 |
| Combined, full default search (0.15–30 d) | 3.18 d not in top 12 | 3.18 d not in top 12 |

Single-sector results are unchanged, except that sector 32's noise is now estimated per orbit. In the combined data, the planet's narrow-window score rises and its depth moves toward the reported ~500 ppm, while the inverted-flux control stays near 7.5–7.8. **The full default search still misses the planet.** Long-period, two-to-four-epoch windows score 19–22 in both the real and the inverted combined data. The current correlated score is a first-pass red-noise correction, not a calibrated false-alarm model. Supplying a stellar density can raise the duration floor and reduce the required period-grid density under its stated geometry assumptions; the cadence-only mode remains available for more inclusive searching. Runs that hit the 250,000-trial cap still cannot guarantee the requested period resolution.

Runtime before → after: combined full search 76.0 → 41.0 s; combined narrow 25.8 → 12.6 s; combined narrow with trend 1.0, 410 → 310 s; serial combined narrow 328.7 → 158.8 s. Single-sector runs stayed within run-to-run variation. Segmentation itself costs milliseconds. The period grid is unchanged, so the speedup comes from how much vetting is triggered; it was not separately profiled. Serial and parallel outputs were byte-identical for the full sector-6 search and the combined narrow search, and the tests compare floating-point bits on mixed-cadence data.

The injection harness (`examples/mixed_cadence_injection.rs`) simulates a 30-minute sector at 250 ppm and a 200-second sector 90 days later at 750 ppm, so both have equal noise per unit time. Each has a 1-day orbit gap. It injects 2-hour box transits at random periods of 1.5–6 days, searches 1–8 days, and counts a recovery when the top periodic dip's period drift over the baseline and its phase offset are both within one duration. The same source runs against both builds.

| Depth | Before: recovered / 12 | After: recovered / 12 | Before → after median score |
| ---: | ---: | ---: | ---: |
| 150 ppm | 0 | 3 | — → 6.30 |
| 200 ppm | 1 | 9 | 6.45 → 6.72 |
| 250 ppm | 2 | 8 | 7.30 → 7.91 |
| 300 ppm | 7 | 11 | 7.49 → 8.67 |
| 400 ppm | 10 | 12 | 8.88 → 10.72 |
| **Total** | **20 / 60** | **43 / 60** | |
| Noise-only nulls with any periodic dip | 0 / 20 | 1 / 20 (score 6.46) | |

Elapsed time was 151 s before and 147 s after for the 80 light curves. Correct weights raise every score, including noise, so the uncalibrated score threshold of 6 now admits one weak null out of 20. Twenty nulls bound the false-alarm rate only loosely (95% interval about 0.1–25%). Calibrating that threshold, for example with signal detection efficiency or larger null sets, is separate work.

## Correlated-noise scoring validation

> Measured with the same earlier defaults as the mixed-cadence section (0.5-day median trend, one search pass, ephemeris deduplication); pass those options to reproduce.

`--scoring correlated` (the default) estimates, for each segment, the correlated-noise factor β at eight log-spaced durations from two cadences to 0.35 days (Pont, Zucker & Queloz 2006). It bins the detrended transit series at each duration with four phase offsets, measures the robust scatter of the bin means divided by their white-noise standard errors, and clamps β to between 1 and 10. Segments with fewer than ten bins borrow the pooled same-cadence estimate, then the whole curve. During the box scan, each window's in-transit variance becomes Σ w·β²/W², using every sample's own segment β interpolated in log duration at that window's duration. Prefix sums of w·β² at each duration knot keep this O(1) per window. Because β ≥ 1, the correlated score never exceeds the white score, so the white score prunes windows first, and the correlated score decides the best window at each period, which windows are discarded, and the ranking before the `--transit-candidates` cap (12 by default). `--scoring white` reproduces the previous statistic exactly. It gave byte-identical candidates to the post-segmentation build in all five TIC 4206066 comparison runs (`scripts/compare_scoring.sh`).

Every periodic dip reports both `score` (the selected statistic) and `white_score`, and correlated runs add `noise_factor` (white ÷ correlated for that window). The CLI prints β per segment, and `--periodogram FILE` writes each trial period's best accepted window. `--min-transit-score` (default 6) and `--transit-candidates` exist for calibration and diagnosis.

**Why the 3.18 d signal is missed** (`scripts/diagnose_missed_signal.py`, combined TIC 4206066 data):

| Statistic | Full oracle | Period oracle | Grid best | Duration-grid loss | Period-grid loss | Uncapped rank |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| white | 15.78 | 15.43 | 15.42 | 0.35 | 0.01 | 32 |
| correlated | 12.77 | 12.25 | 12.23 | 0.52 | 0.01 | 35 |

The oracles are dense searches near the signal: 20,000 periods, with 64 or 10 durations. Grid resolution costs about 2–4% of the score, and the capped period grid costs almost nothing, so **the miss is a ranking loss**. Correlated scoring does not fix it. The planet's 28 transits include the noisiest segment (sector 32's second orbit, β ≈ 1.6–1.8 at 1–2 h), so its noise factor is 1.26. The stronger peaks are 2–3-event windows at 16–30 days whose events fall in quiet segments, with factors of 1.07–1.24. They come from a few localized excursions, which a robust β is deliberately insensitive to. Many are near-duplicates of the same two-event alias (for example, 28.82–28.85 d), because a two-event period is so loosely constrained that drift-based deduplication keeps each copy. An event-structure test and deduplication by shared events are the next ranking fixes. Physical duration limits and period-grid changes are a separate follow-up.

**Matched false-alarm comparison** (`examples/correlated_noise_calibration`, 1,706 s; results in `benchmarks/correlated-noise-calibration/`):

- **Simulations.** A 30-minute sector (250 ppm per sample) and a 10-minute sector 90 days later (433 ppm per sample), each with two 13.2-day orbits. In the correlated model, both get an Ornstein–Uhlenbeck process (60 and 160 ppm; 0.08 and 0.12 d), scattered-light ramps at orbit starts, and about 1.5 Gaussian dip-like glitches per orbit, none of which the quoted errors include.
- **Real data.** The 48 sector-32 TESS-SPOC curves in `data/calibration/`.
- **Injections.** Trapezoid transits at fixed white-noise SNR. A detection requires the top periodic dip to match the injection (or any listed dip, for two-planet and real-curve families) with a score strictly above the null threshold.

| Null set (white / correlated statistic) | Curves | T at 5% FAR | T at 1% FAR | Curves ≥ 6.0 |
| --- | ---: | ---: | ---: | ---: |
| White noise, P 1–10 d | 300 | 5.90 / 5.88 | 6.18 / 6.14 | 6 / 6 |
| White noise, P 20–50 d | 150 | 5.80 / 5.78 | 6.13 / 6.13 | 5 / 5 |
| Correlated noise, P 1–10 d | 300 | 6.87 / 6.04 | 7.30 / 6.39 | **195 / 17** |
| Correlated noise, P 20–50 d | 150 | 6.95 / 5.95 | 7.28 / 6.28 | **68 / 6** |
| Inverted real curves, P 1–10 d | 48 | 119.4 / 22.9 | 525.6 / 56.2 | 32 / 26 |

| Recovered (white → correlated) | Injections | 5% FAR | 1% FAR |
| --- | ---: | ---: | ---: |
| White noise: 1 / 2 / 4 / 6 h | 36 each | 32→32 / 33→32 / 26→25 / 16→16 | 30→30 / 29→28 / 24→24 / 14→14 |
| White noise: two-event; two-planet inner / outer | 45; 40 / 40 | 37→37; 40→40 / 40→40 | 35→34; 40→40 / 40→40 |
| Correlated noise: 1 / 2 / 4 / 6 h | 36 each | 24→25 / 21→21 / 18→16 / 13→11 | 22→24 / 16→15 / 16→14 / 11→6 |
| Correlated noise: two-event; two-planet inner / outer | 45; 40 / 40 | 24→22; 40→40 / 36→34 | 24→22; 39→39 / 32→26 |
| Real curves, per-star inverted null: 2 h / 4 h | 96 each | 47→53 / 41→42 (one threshold per star) | |

What this shows:

- **Calibration.** With the same fixed threshold of 6, the white statistic flags 65% of correlated-noise null curves and the correlated statistic flags 5.7%, while white-noise nulls are unchanged (2% either way). The correlated score's meaning therefore varies much less with the noise type.
- **Power at matched FAR did not improve in these simulations.** Short transits gain slightly. Long transits, two-event planets and the outer planet lose some, most clearly 6-hour transits at 1% FAR (11 → 6). With 36–45 injections per row, differences of two or three are within binomial noise (about ±3). Penalizing long windows lowers real long transits as much as long noise windows, and the dominant simulated false alarms are glitches, which β does not target.
- **Real curves.** On real curves with per-star thresholds, 2-hour recovery improves (47 → 53 of 96). Pooled real thresholds are dominated by inverted variable stars and eclipsing binaries, so they are not a usable false-alarm calibration. Survey use needs variability screening or per-star nulls.
- **Determinism.** 24 serial-versus-parallel spot checks in the harness matched exactly, and `tests/correlated_noise.rs` compares floating-point bits for candidates, periodograms and β profiles.

Runtime, measured before other concurrent changes to the search loop: correlated scoring took 0.7–1.5× the white time on the TIC 4206066 runs (full combined 37.2 vs 44.4 s; narrow combined 17.9 vs 12.4 s), and 1.9× in the calibration harness (977 vs 522 s). Its extra work is the per-knot histograms, and fewer windows passing the threshold offsets part of that.

## Mac Studio performance

Measured locally on a Mac Studio with M2 Ultra, 24 CPU cores (16 performance / 8 efficiency), 60 GPU cores, and 128 GB memory. The table gives median end-to-end elapsed time across three runs of the included 1,350-observation synthetic fixture, with period bounds 0.5–10 days and `periods=1200`. All compared CSV outputs were byte-identical to the previous serial implementation. The 20 tests also pass for portable and native builds, including a comparison of candidate floating-point bits between serial and parallel scans.

| CPU workers | Portable release | Native M2 release |
| --- | ---: | ---: |
| 1 | 2.819 s | 2.829 s |
| 8 | 0.426 s | 0.412 s |
| 16 | 0.273 s | 0.273 s |
| 24 | 0.240 s | 0.230 s |
| Automatic | 0.249 s | 0.253 s |

The portable 24-worker build was about 11.7× faster than serial for this fixture. Native CPU targeting did not show a consistent advantage beyond run variability, so the normal portable release remains a useful default. Full measurements are in `benchmarks/mac-studio-m2-ultra.json`; larger or gapped multi-sector data need their own timing and recovery checks. The benchmark demonstrates CPU scaling, not improved scientific completeness.

```sh
cargo build --release
/usr/bin/time -p ./target/release/astrosearch \
  --input benchmarks/fixtures/transit-1350.csv \
  --min-period 0.5 --max-period 10 --periods 1200 --threads 24 \
  --output /tmp/astrosearch-benchmark.csv
```

For an optional build tuned to the current host, keep the output separate:

```sh
RUSTFLAGS='-C target-cpu=native' CARGO_TARGET_DIR=target/native cargo build --release
./target/native/release/astrosearch --input target.csv --threads 24
```

`target-cpu=native` specializes code for the build host, so use the normal release when portability is required ([Rust compiler documentation](https://doc.rust-lang.org/rustc/codegen-options/index.html#target-cpu)). The GPU is currently unused. A future Metal path could batch a coarse search in 32-bit arithmetic and refit shortlisted candidates on the CPU in 64-bit arithmetic, but would need explicit precision, recovery, and timing validation: Metal does not provide a native `double` type ([Metal precision limitation documented by PyTorch](https://docs.pytorch.org/docs/main/notes/mps.html#double-precision-types)). The CPU's 64-bit reference path should remain available for validation.

## Research basis and recommended next steps

### 2026-10 search-grid and runtime update

The transit search now uses a period grid spaced to limit timing drift to half the shortest searched duration. With no stellar-density prior, it retains the previous uniform log-period grid; supplying `stellar_density_g_cm3` redistributes the available trials using a conservative minimum-duration model (impact parameter up to 0.95 and eccentricity up to 0.9). The model is a search prior, not a physical guarantee. If the required count exceeds the configured cap, the search warns and may miss narrow peaks. The per-period repeat-event check now uses sorted per-cycle phase ranges and prefix sums, instead of rescanning every observation for each improved trial window.

On the local mixed-cadence TIC 4206066 combined curve (34,511 valid observations, 2,577.55 days), the preserved pre-change release took 34.50 s with 24 workers and the updated release took 28.96 s in one end-to-end run (about 16% faster). Candidate CSVs were byte-identical on this same input. This is one machine/run comparison, not a general speedup guarantee. The search still hit its 250,000 period-trial cap. A one-injection smoke calibration on Sector 32 recovered 0/1 injected 200-ppm signals in both the old and new builds; that sample is far too small to estimate completeness, but confirms the indexing change did not change that outcome. Run the validation tool with many randomized injections per depth and target before using score thresholds for survey claims.

- The Rust detector remains a fast box-template screen and now masks each strongest signal before repeating the search. An optional iterative TLS search uses limb-darkened transit templates; TLS's white-noise FAP is not a calibrated false-alarm probability. The official [TLS interface documentation](https://transitleastsquares.readthedocs.io/en/stable/Python%20interface.html) describes the template search and its inputs. Compare both engines with injection/recovery across quiet, active, and gapped curves before setting survey thresholds ([Kovács, Zucker & Mazeh 2002](https://arxiv.org/abs/astro-ph/0206099); [Hippke & Heller 2019](https://arxiv.org/abs/1901.02015)).
- The default trend is now a sliding Tukey-biweight estimate, and candidate depths are remeasured after masking transit windows. This is a local implementation inspired by Wōtan's benchmarked filter family, not a direct Wōtan port; compare recovery against the median option for active and gapped stars ([Hippke et al. 2019](https://arxiv.org/abs/1906.00966)).
- Odd/even depths are one false-positive diagnostic, but complete TESS vetting also uses secondary-eclipse/shape tests and pixel-level centroid/difference-image information. The official SPOC pipeline produces calibrated pixels, light curves, transit searches, and validation products; this tool currently consumes only CSV light curves ([NASA HEASARC TESS pipeline documentation](https://heasarc.gsfc.nasa.gov/docs/tess/documentation.html); [LEO-Vetter methods paper](https://doi.org/10.3847/1538-3881/ae070a)).
- Before expanding to a large target sample, measure completeness and false alarms with injection/recovery and time-correlated null controls. Do not turn score thresholds into discovery claims or planet probabilities without those calibrations.

# Historical JEV vs local target-ranking evaluation

This report documents a completed comparison only. JEV has since been removed from the active survey runner; new searches use the local deterministic ranking.

## Design

- Source: MAST TIC filtered cone query, 5° radius centered at RA 80.85984096°, Dec −1.32727827°. MAST returned 229 targets satisfying Teff 2700–5200 K, distance 0.1–100 pc, stellar radius 0.05–0.9 R☉, and Tmag 4–14. This uses MAST’s documented filtered TIC position service ([MAST API](https://mast.stsci.edu/api/v0/_services.html)).
- Both ranking methods received the same 229 target cohort. Jev received TIC ID, Teff, distance, radius, and Tmag. The deterministic local comparator ranked by the approximate photon-limited Earth-transit S/N proxy `10^(-0.2 × (Tmag − 10)) / radius_solar²`. This is a transparent baseline, not a fitted TESS noise model.
- Each method selected its top 10. Both used QLP only, the same 1095-day minimum baseline, at least one sector ≥90, 1–10 day periods, 12,000 requested trials, and 24 Rust threads. SPOC had no matching recent FFI products in the previous pilot, so it was not added to this batch.
- The top-10 sets overlap on 6 TICs. The shared targets were downloaded/searched once and reused in both method comparisons; each method has 4 unique selections.

## Selection and usable-data yield

| method | selected | full multi-year QLP searches | coverage/data outcome |
|---|---:|---:|---|
| Jev | 10 | 4 | 4 below 1095 days; 1 without sector metadata; 1 without usable products |
| local S/N proxy | 10 | 6 | 1 below 1095 days; 1 without sector metadata; 2 without usable products |

Jev therefore yielded 4 searchable targets out of 10; the local proxy yielded 6 out of 10. On this cohort, Jev’s metadata-only ranking did not account for whether QLP had usable multi-year coverage.

## Injection/recovery

We injected 1.0 R⊕ and 1.5 R⊕ box transits into each of the 7 distinct selected targets with searchable QLP curves. The two test periods were 2.37 and 6.43 days, with fixed phase fractions 0.17 and 0.61; transit depth was scaled to each TIC stellar radius. We reran the same Rust search settings. A recovery means a listed periodic-dip candidate matched both the injected period and epoch within one estimated transit duration. Each unmodified candidate list served as a same-period/phase null check.

| method | injected planet radius | recovered | same-window matches in unmodified curves |
|---|---:|---:|---:|
| jev | 1.0 R⊕ | 2/4 (50.0%) | 0/4 |
| jev | 1.5 R⊕ | 2/4 (50.0%) | 0/4 |
| local | 1.0 R⊕ | 3/6 (50.0%) | 0/6 |
| local | 1.5 R⊕ | 4/6 (66.7%) | 0/6 |

No same-window matches occurred in the 14 unmodified-curve controls. Recovery was 50% for Jev versus 50% for local at 1.0 R⊕, and 50% versus 66.7% at 1.5 R⊕. The small sample produces wide uncertainty; it does not establish a statistically reliable difference or survey-wide completeness. Recoveries may be listed below the overall candidate rank 1, so they are search triggers rather than validated planets.

## Cost and runtime

- Jev made 5 API requests for the 229-target ranking, using 59,053 input tokens and 3,803 output tokens. At the published $0.042 per million input tokens and $0 output-token rate, estimated inference cost was **$0.002480** (about 0.25 cents).
- The Jev arm took 218.5 seconds for model ranking plus 10 selected targets. The local arm took 232.0 seconds to process its 4 unique targets after reusing the Jev/local overlap. These elapsed times are not a full apples-to-apples batch-speed comparison. Jev does not speed up the Rust search; it adds a small ranking cost.

## Readout

On this 229-target cohort, Jev did not improve the target yield over the deterministic baseline. It selected fewer targets with usable multi-year QLP data (4 vs 6); the injection sample showed no improvement at 1 R⊕ and a modest, inconclusive local advantage at 1.5 R⊕. Given the low cost, Jev is affordable as an optional ranking, but these results do not justify paying for it to improve search performance. The better next improvement is to include archive coverage and measured light-curve noise before using the target cap.

## Artifacts

- `decision-eval-targets.csv`: MAST TIC cohort.
- `decision-eval-v2-rankings.csv`: complete method ranks and Jev scores.
- `decision-eval-v2-selection-results.csv`: ranks joined to QLP coverage outcomes.
- `decision-eval-v2-jev/` and `decision-eval-v2-local/`: selected target products, manifests, candidate CSVs, and run ledgers.
- `decision-eval-v2-injections/injection_results.csv`: per-injection outcomes; `injection_summary.json`: grouped totals.
- `decision-eval-v2-jev/survey-ledger.sqlite` and `decision-eval-v2-local/survey-ledger.sqlite`: per-arm cross-run checkpoints.

**Limits:** this is one 5° sky cone, one TIC catalog snapshot, two deterministic injection periods/phases, and only the targets selected by these two methods. TIC catalog inputs, QLP pipeline products, and an approximate box-transit injection do not replace blinded injection/recovery across many fields, randomized phases, period/depth grids, null curves, or independent-pipeline validation.

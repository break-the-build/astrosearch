# Signal diagnostics: depth provenance, event patterns and variability conflicts

Branch `diag/variability-depth-20261011`, based on `main` at `6d9f428`. TIC 366989133 (WT 682) served as an adversarial case, not as ground truth: its nature is **unresolved**. Catalogue checks failed (DNS/API errors), no target-pixel localization exists, and nothing here classifies it as a planet or as a known false positive.

## What the WT 682 outputs showed

These come from the Profile J five-pass SPOC search in the live survey (`profile-j-run-20261010/.../TIC366989133`, read-only; reproduced from a copy).

- The top `periodic_variability` signal has period 0.37505136 d, amplitude 3,577 ppm and score 199.5.
- Three transit-like rows reached the 25.37 detection tier: 0.18772020 d, 0.18769245 d and 0.37544072 d.
- The two near-identical 0.1877 d rows report depths of 1.89 ppm and 2,694.9 ppm.
- TLS returned the same 0.37544 d feature on three iterations and labelled the repeats "integer" aliases.

## Root causes, from the code

1. **The depth mismatch is an estimator switch, not a unit error.**
   - `src/refine.rs::remeasure_depths` replaces the box-fit depth with a local depth: the median of samples 1.5–3 durations away minus the median in transit, after transit-masked detrending.
   - It does this only when the local value is positive, and it overwrote the box depth without keeping it.
   - On a smoothly varying star the local value scatters around zero. Row 21's local depth was +1.9 ppm, so it replaced the 2,837.1 ppm box depth. Row 22's was −115.7 ppm, so its 2,694.9 ppm box depth was kept.
   - `depth_ppm` therefore mixed two estimators depending on the sign of a noise-level number.
2. **The half-period detection-tier rows are created by multi-pass masking.**
   - With one pass, the same search returns no 0.1877 d rows. All 12 transit candidates are near 0.375 d, with scores of 16.1–18.1, all below the tier.
   - With five passes, the 0.1877 d rows come from pass 5, after 55.2% of samples had been masked. The 0.37544 d detection-tier row comes from pass 3, after 28.7% had been masked.
   - Masking the modulation minima raises the out-of-transit mean of what remains. Row 21's events then measure an average of 2,815 ppm against that baseline, while the same events average 45 ppm in the unmasked light curve (event-by-event correlation 0.92).
3. **Rows describing one feature were not linked.**
   - The Rust event grouping merges rows that share at least 80% of the smaller event set. Rows 21 and 22 share 79% (277 events), so they stayed separate.
   - Most of the 0.1877 d rows share 44–80% of their events with each other.
4. **TLS repeats and units.**
   - TLS fits a 0.0019 d duration to a dip about 0.1 d wide, so masking the fitted window leaves the feature in place and the next iteration refits it at a shifted epoch. Masking twice the duration did not change this on WT 682, so the mask is unchanged.
   - A ratio of 1 was labelled "integer".
   - The committed script also wrote TLS's bottom-flux level as a depth (996,135 "ppm" instead of 3,864.6 ppm). Codex's uncommitted working copy has the same one-line fix.

## Changes

None of these changes alters a score, period, epoch, depth value, tier, threshold or the candidate list. On WT 682 the five-pass output from this branch matches the survey output in every cell of all 130 rows except `note`.

- **`src/refine.rs`:** `note` keeps the replaced box depth: `depth remeasured after transit-masked detrending (box depth N ppm)`.
  - A local depth below 25% of the box depth adds a `depth inconsistency` warning.
  - A non-positive local depth is stated (`transit-masked local depth N ppm is not positive; box depth kept`).
  - The warning only works one way. On synthetic sinusoids, local/box ratios ranged from negative to 0.95 across rows, so the absence of the warning is not evidence of a sharp transit. An injected box transit gave 1.02.
- **`src/transit.rs`:** candidates from pass 2 onward record `found in search pass k with x% of samples masked by earlier passes`.
- **`scripts/vetting/signal_diagnostics.py` (new) and `scripts/signal_diagnostics.py` (new CLI, `astrosearch-signal-diagnostics`):** report-only diagnostics, written to new JSON and Markdown files.
  - Event patterns: links at ≥50% shared observed events. The match tolerance is max(half the summed durations, 2% of the shorter period). Each link is labelled with its period relation.
  - Depth check: every depth with its estimator, plus the depth on the unmasked light curve at the same ephemeris. Disagreements are flagged; no value is chosen.
  - Variability conflict: for periods at 1/2, 1 or 2 times a reported modulation period (0.2% tolerance), two harmonics are fitted per segment to out-of-transit data and the dip is measured before and after. Below 30% retained is reported as a competing explanation; above 70% as not explained by the modulation.
  - Interleaved pairs: events at a consistent phase of one period. A phase within 0.05 of 0.5 is reported as primary/secondary-like; other offsets only beyond 0.2 of the period. Each row's parity coverage is also reported.
- **`scripts/tls_search.py`:**
  - Correct depth units (`depth_ppm_from_tls`).
  - A new `observed_transit_times` column, and `shared_events_with_earlier`, which lists links at ≥50% shared events.
  - The `same-period` label.
  - Search results are unchanged; the mask stays at 1× the duration.
- **`scripts/evaluate_signal_diagnostics.py` (new):** the flag-rate evaluation described below.
- **Tests:** `tests/test_signal_diagnostics.py` (unit tests plus synthetic Profile J five-pass cases: transit, sinusoid, eclipsing binary with unequal eclipses, tenfold-diluted blend, noise), `tests/test_tls_dedup.py`, and one Rust regression test.

## WT 682 after the changes

From `astrosearch-signal-diagnostics` on the reproduced five-pass output plus the TLS rows. The report is in `data/diagnostics/signal-diagnostics-20261011/wt682-diagnostics.{json,md}`; its input paths point to a scratch copy of the survey files.

- **Event patterns:** 11 astrosearch rows form one event pattern (0.1877 d and 0.375 d rows, linked as same-period and 2:1).
  - The local depth is at most 23% of the box depth in every row, or negative.
  - Ten rows have unmasked depths of −18.6 to 438.6 ppm against box depths of 1,948–2,947 ppm.
- **The 0.37544 d detection-tier row:** a separate pattern sharing no events with the 0.1877 d rows. Its local depth is 18% of its box depth.
- **TLS iterations** (the survey's three-iteration file): they form their own pattern. They are linked as same-period with 61–99% shared events, and their depths agree (3,854–3,891 ppm).
- **TLS rerun on this branch:** four iterations with arguments identical to `main`'s script (`--min-period 0.15 --max-period 1.0 --max-planets 4`, stellar radius 0.208428, mass 0.177393).
  - Every search column is identical to `main`.
  - `depth_ppm` changes from about 996,135 to 3,838.8–3,890.6 ppm.
  - Iterations 2–4 are labelled `same-period`. They share 178/291, 186/237 and 242/247 events with an earlier iteration.
- **Modulation test:** all 15 transit-like rows have related modulation at 0.37505 d.
  - 14 are reported as explained by it.
  - One (row 32) had no positive unmasked dip to test.
  - For the strongest, removing two harmonics leaves 1%: astrosearch 0.37544 d goes from 2,977.9 to 35.7 ppm, and TLS iteration 1 from 3,544.2 to 40.5 ppm.
- **Reading:** the transit-like detections are explained by the star's periodic modulation, as measured on the unmasked data. This is a competing explanation, not a classification. Spotted rotation, ellipsoidal or reflection variation, or a contact-binary-like light curve all remain possible, and the signal's origin on the sky is unverified.

## Synthetic regression cases (Profile J, five passes)

Two 13-day segments at 10-minute cadence, with 400 ppm white noise.

- **Transit** (2.3 d, 2,500 ppm): recovered in pass 1, depths consistent, not explained by modulation.
- **Pure sinusoid** (0.6 d, 3,000 ppm): every tested related dip retains less than 30% after modulation removal.
- **Eclipsing binary** (1.8 d, 9,000 and 3,000 ppm eclipses):
  - Pass 1 finds the 1.8 d primary. Pass 2 reports the secondary as a 0.90 d row whose events all fall on odd cycles, so its odd/even statistic is blank.
  - The interleaved-pair check links them at phase 0.502, with box depths of 8,948.6 and 2,964.7 ppm (ratio 3.02).
- **Tenfold-diluted blend:** the same structure, with depths of 848.6 and 333.5 ppm.
- **Noise:** the pipeline runs cleanly.

## Flag rates on real light curves

Measured with `scripts/evaluate_signal_diagnostics.py`, using this branch's binary with Profile J and five passes (`--search-profile fast --transit-search-passes 5`, other options at their defaults, no per-star density). Cases are the known-answer **calibration split only**; the blind split was not read. The real eclipsing binaries are the 132 curves of the hard-negative survey (`data/benchmarks/hard-negatives/survey-final`). Results are in `data/diagnostics/signal-diagnostics-20261011/` (`summary.json`, `cases.csv`).

```sh
uv run python scripts/evaluate_signal_diagnostics.py \
  --known-answer data/benchmarks/generated-20261009/truth-manifest.csv \
  --eb-survey data/benchmarks/hard-negatives/survey-final \
  --binary target/release/astrosearch --output <dir> --workers 14
```

For planets the flags are read on the candidate that matches the catalogue ephemeris. For controls and binaries they are read on the highest-jackknife transit candidate. The selections differ, so the columns are not a matched comparison. Percentages are of evaluated cases; brackets give 95% Wilson intervals.

| Flag | Planets (matched), n=175 of 200 | Inverted controls (strongest), n=178 of 200 | Real eclipsing binaries (strongest), n=131 of 132 |
|---|---|---|---|
| Reaches the 25.37 detection tier | 74 (42%) | 7 (3.9%) | 124 (95%) |
| Explained by related modulation | 0 [0–2.1%] | 14 (7.9%) | 1 (0.8%) |
| Interleaved at phase ~0.5 (secondary-like) | 0 [0–2.1%] | 1 (0.6%) | 30 (23%) [17–31%] |
| Events on one cycle parity only | 0 [0–2.1%] | 10 (5.6%) | 1 (0.8%) |
| Depth depends on masking (unmasked < 25% of box) | 7 (4.0%); 0 of 74 at tier | 23 (13%) | 0 |
| Depth inconsistency (any depth-check issue) | 8 (4.6%); 6 of 74 at tier | 38 (21%) | 55 (42%) |

Excluded cases: 25 planets had no candidate matching the catalogue ephemeris, 22 controls and 1 binary had no transit candidate.

What this does and does not show:

- **Planets:** the modulation, interleaving and parity flags never fired on a matched real-planet candidate in this sample (upper 95% bound 2.1%). The masking flag fired on 7 matched planets, none of them at the detection tier. The depth-inconsistency flag fired on 8, including 6 at the tier, so it is a prompt to read the depths, not evidence against a candidate.
- **Eclipsing binaries:** 95% of them reach the detection tier, so the tier does not reject them. The interleaving check marks 23% (30 of 124 at the tier) with a secondary-like pair, independently of the odd/even statistic.
- **Not a detection result:** nothing here is a detection-performance or false-alarm claim. The diagnostics are report-only: they remove no candidate and change no threshold, so recall and the false-alarm rate at the tier are unchanged by construction. A rule that acts on these flags would need its own matched-false-alarm evaluation on held-out data.
- **No search behaviour changed:** on 20 benchmark cases (10 planets, 10 controls; 766 rows), this branch's binary and `main`'s gave identical cells in every column except `note`, and every original note fragment was preserved.

## Limitations and follow-up

- **Multi-pass masking:** on strongly variable stars it can create detection-tier rows (WT 682). The notes now expose it, but the search still does it, and the five-pass Profile J tier was calibrated with this behaviour. Possible fixes include re-estimating the baseline excluding only the current pass, stopping passes when a related modulation dominates, or prewhitening the strongest modulation before searching. Each changes the search and needs a matched-false-alarm evaluation with eclipsing-binary, blend and systematic negatives.
- **The Rust 80% grouping rule is unchanged.** Lowering it would change the candidate list the tier was calibrated on. The Python layer reports links at 50% instead.
- **The modulation test only runs for rows near a reported `periodic_variability` period.** It uses a fixed two-harmonic model per segment, so non-sinusoidal or evolving variability may be under-fitted, and very broad transits could be partly absorbed.
- **The interleaved-pair check needs both eclipses to be reported as rows.** A secondary below the screening floor will not appear.
- **Pre-existing on `main` (`6d9f428`):** six Rust tests fail, all of which assume the old standard default configuration:
  - `transit::tests::duration_dependent_beta_changes_which_window_is_selected`;
  - `sectors_report_their_own_cadence_and_noise`;
  - `two_stage_search_recovers_the_full_search_period_with_fewer_trials`;
  - `two_transits_can_have_period_longer_than_half_baseline`;
  - `max_period_trials_caps_the_grid_and_warns`;
  - `fast_folding_coarse_stage_recovers_the_full_search_result`.

  `src/main.rs` is also not `rustfmt`-clean. This branch does not touch either.
- **Not done:** wiring the diagnostics into `tess_survey.py` and the CTOI gate. Codex has large uncommitted changes to `tess_survey.py`, so this is left as follow-up to avoid conflicts.
- **WT 682 still needs:** catalogue/TCE/variable-star queries once connectivity returns, Gaia RUWE and neighbours, and target-pixel localization.

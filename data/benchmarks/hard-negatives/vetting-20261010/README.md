# Light-curve vetting against planets and false positives (2026-10-10)

Produced by `scripts/evaluate_vetting_negatives.py` (blind split, seed 20261010, `--planets 100 --blends 60`). `vetting_negatives.json` is the summary; `vetting_negatives.csv` has one row per case (388).

**Search and threshold.** All candidates come from the **standard** search profile, and the detection tier is the standard-profile jackknife threshold **23.30** (`data/calibration/detection_threshold.json`). These results do not apply unchanged to Profile J, whose recalibrated threshold is different; rerun with `--detection-calibration` for J.

**Checks.** Light-curve vetting = secondary-eclipse scan, odd/even depth test, and grazing shape from the physical fit. It does not include the gate's pixel-centroid, TRICERATOPS, or empirical-FAP stages.

| Class | Source | Pass light-curve vetting | Pass vetting and detection tier |
|---|---|---|---|
| Confirmed planets | known-answer benchmark | 62 / 79 vetted (78%; Wilson 95% 68–86%) | 29 / 79 (37%) |
| Real eclipsing binaries | TESS EB catalogue (Prša et al. 2022), 128 of 136 searched | 14 / 128 (11%; 7–18%) | 12 / 128 (9%; 5–16%) |
| Systematics | inverted light curves | 61 / 99 (62%) | 3 / 99 (3%) |
| Blends | synthetic diluted EBs | 27 / 60 (45%; 33–58%) | not computed |

Rejections:
- **Planets:** odd/even 14, grazing 2, secondary 1.
- **Eclipsing binaries:** odd/even 73, grazing 69, secondary 58 (one case can fail several checks).
- **Blends:** secondary 21, odd/even 12, grazing 4.

**Reading.** Light-curve vetting removes about 89% of real eclipsing binaries, but about 9% still pass both it and the detection tier. The tier barely filters eclipsing binaries, because their signals are real. About half of the synthetic blends pass. The pixel and statistical-validation stages are therefore required before any candidate is treated as planet-like. The odd/even test also rejects 18% of real planets, which suggests its threshold should account for correlated noise.

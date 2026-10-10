# TIC 4206066 TESS light curves

Real TESS photometry of TIC 4206066 (StKM 1-561), used as a known-answer check.
Rabtsevich (2026) reports a ~500 ppm transit candidate at P = 3.18 d (23 transits)
and a tentative 11.13 d signal ([write-up](https://doi.org/10.5281/zenodo.22967456),
[data and code](https://doi.org/10.5281/zenodo.22999133)). Neither is confirmed.

| File | Sector | Product | Cadence | Rows | Median quoted error |
| --- | ---: | --- | ---: | ---: | ---: |
| `s06_TESS-SPOC_1800.csv` | 6 (2018–19) | TESS-SPOC FFI | 1800 s | 975 | 217 ppm |
| `s32_SPOC_120.csv` | 32 (2020) | SPOC target pixel | 120 s | 17,094 | 840 ppm |
| `s98_QLP_200.csv` | 98 (2025) | QLP FFI | 200 s | 16,442 | 831 ppm |
| `combined.csv` | all three | time-sorted concatenation | mixed | 34,511 | — |

Columns are `time` (BTJD), `flux`, and `flux_err`. Each sector was downloaded from
MAST with lightkurve 2.6.0 (numpy 2.5.3, astropy 8.0.1) using
`quality_bitmask="hard"`, NaNs removed, and flux normalized to its own median, so
sectors share a unit baseline. No other detrending was applied.

```sh
uv venv .venv && uv pip install --python .venv/bin/python lightkurve==2.6.0
.venv/bin/python -I scripts/export_tess.py data/tic4206066 "TIC 4206066" \
  6:TESS-SPOC:1800 32:SPOC:120 98:QLP:200
```

The export also writes a `cache/` directory of downloaded FITS files, which is not
kept here. SHA-256 checksums of the committed CSVs:

```
c86f5c6a33a5859d718a576ef684caad7706259400c11753d350bf02ac5f1533  s06_TESS-SPOC_1800.csv
0a372977cba9a3b28e3ec85dce0199ee1c033b6a612c64c698089e8cc07c20ba  s32_SPOC_120.csv
8c753b6663fd71f851a15ebdc955622ae8de86b9da3fcf4c78790da2099c0ff0  s98_QLP_200.csv
f62a6163690dbe899e0deabce8dc8bb1ac11a9c033fbfe10319f623c8e3dae7b  combined.csv
```

`scripts/review_runs.sh` reproduces every search, null control, and serial/parallel
comparison reported in the top-level README.

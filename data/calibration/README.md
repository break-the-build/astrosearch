# Calibration light curves

`s32_tess_spoc_600/` holds 48 real TESS-SPOC FFI light curves from sector 32
(600-second cadence, 164,073 rows in total): the 48 stars nearest TIC 4206066,
within 1,308 arcsec, that have such a product. They share the target's camera
and its scattered-light and pointing systematics. Several are variable stars or
eclipsing binaries; none was removed, so the set includes real astrophysical
false-alarm sources.

`examples/correlated_noise_calibration` uses them two ways:

- **Inverted nulls.** Flux is reflected about its median, so real dips become
  brightenings and any periodic dip found afterwards comes from noise,
  systematics, or inverted variability.
- **Injections.** 2- and 4-hour trapezoid transits are injected into the
  original curves at fixed white-noise SNR.

Each file has `time` (BTJD), `flux` (normalized to its median), and `flux_err`,
written with `quality_bitmask="hard"` and NaNs removed. `manifest.csv` lists
each TIC ID, row count, angular distance, and SHA-256 checksum
(`manifest.csv` itself: `5c5d2f8ca6dfcb52a21bdee4622b46da5ca2c6ef841a2b32f39c355a9ade3eb7`).

```sh
uv venv .venv && uv pip install --python .venv/bin/python lightkurve==2.6.0
.venv/bin/python -I scripts/fetch_calibration_curves.py \
  data/calibration/s32_tess_spoc_600 "TIC 4206066" 32 600 3600 48
```

MAST search results can change as products are added, so a new download may
select different stars; the committed files and checksums are the reference set.

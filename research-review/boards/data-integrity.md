# Board 4 — Astronomical data integrity and calibration

**Sequential simulated review; not an independent external audit.**

## Finding

**Rating: 3.5/5 — strong light-curve provenance for this release; source-local
validation remains incomplete.**

**VERIFIED:** benchmark provenance records the NASA Exoplanet Archive query,
catalog snapshot checksum, seeded split, product identifiers, catalog
ephemerides, compact-curve hashes, Lightkurve version, and exclusions. The
auditor checked the 600 actual run inputs against the truth metadata and every
candidate CSV against its sidecar. All checksums matched. The builder excludes
uncertain/TTV-flagged ephemerides and verifies observed events before inclusion.

The sample contains 3,824 QLP products in Sectors 1–104. The MAST QLP page says
QLP uses TESS-Gaia Light Curve methodology starting Sector 94; 276 sample product
sectors fall at 94–104. Therefore QLP is not one homogeneous extraction
pipeline across the release. This must be recorded when comparing search
versions.

**UNKNOWN:** whether every light curve's source time system, flux column,
normalization, quality mask, and archived pipeline revision were preserved
through all transformations. The benchmark's compact CSV is traceable to
products, but it does not retain pixels or support source localization.

## Recommendations

- Retain product filenames, sector, author, processing version, selected flux
  column, quality mask, original time convention, and normalization per sector.
- Add round-trip checks for BJD/BTJD, exposure time, units, and cadence
  transitions on the product families actually downloaded.
- Stratify Sectors 94+ separately or normalize extraction-family changes using
  common targets when comparing algorithm performance.
- Keep the pixel-level provenance chain for candidate vetting; a light curve
  alone cannot establish that the transit is on the target.

Official references: [MAST QLP data products](https://archive.stsci.edu/hlsp/qlp),
[MAST TESS-SPOC products](https://archive.stsci.edu/hlsp/tess-spoc), and
[NASA TESS pipeline documentation](https://heasarc.gsfc.nasa.gov/docs/tess/data_pipeline.html).

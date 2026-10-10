# Board 1 — Scientific review

**Simulated assessment; not written by or endorsed by the named researchers.**
Evidence labels follow the program brief.

## Finding

**Rating: 3/5 — adequate as a photometric candidate-triage prototype.**

**VERIFIED:** Search input is a time, positive flux, and optional error series;
the Rust program uses segment-aware detrending, correlated-noise weighting,
box-shaped periodic dips, sinusoidal variability, and isolated events. The
public TESS QLP/SPOC front end and downstream vetting tools are real code, but
the known-answer evaluation tested only QLP one-dimensional curves and one
transit-search pass.

**SUPPORTED:** A box model is a valid computational screening model for
periodic transit-like events. TLS research supports limb-darkened templates as
a useful small-planet detector and independent check; Wōtan benchmarks support
time-windowed robust biweight detrending. These papers support trying those
methods, not assuming the local versions reproduce their reported performance.
The [NASA SPOC pipeline](https://heasarc.gsfc.nasa.gov/docs/tess/data_pipeline.html)
also searches calibrated TESS data and produces TCE/DV products, so the local
tool should be positioned as an independent search and triage path.

The benchmark recovered 87/200 known blind planets at its calibration-selected
jackknife threshold (43.5%, Wilson 95% CI 36.8–50.4%). Recovery ranged from
15.4% for shallow-depth positives (4/26) to 75.4% for deep positives (46/61).
These figures do not imply completeness for the TESS population. The 583/600
period-grid cap warnings make the search resolution a direct scientific
limitation.

## Recommendations

1. Keep the result label “transit-like candidate” through all automated stages.
   Do not call a signal a planet or discovery without independent validation.
2. Treat the period grid cap as a sensitivity defect: make the grid complete
   for declared bounds or emit a candidate-quality flag and measure recovery
   against an uncapped oracle/injection suite.
3. Compare Rust box search, TLS, and SPOC/QLP TCEs on one frozen target set at
   matched empirical false-alarm rates. Preserve the calibration split.
4. Group results by sector/product cadence and QLP processing family. Sectors
   94+ use changed QLP/TGLC methodology per MAST documentation.
5. Expand negative controls to real eclipsing binaries, known false alarms,
   background blends, and time/block-shift nulls. Inverted light curves alone
   cannot establish astrophysical reliability.

**UNKNOWN:** recovery under TTVs, eccentric/grazing transits, transit timing
variations, variable-star activity, sector boundary effects, and true
multi-planet overlap at scale. The benchmark has no independent planet
confirmation experiment.

Sources: [TLS](https://arxiv.org/abs/1901.02015),
[Wōtan](https://arxiv.org/abs/1906.00966),
[SPOC overview](https://heasarc.gsfc.nasa.gov/docs/tess/data_pipeline.html).

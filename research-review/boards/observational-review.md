# Board 7 — Observational follow-up and discovery confirmation

**Sequential simulated review; not an independent external audit.**

## Finding

**Rating: 2/5 — workflow components exist, but no new candidate has cleared an
observational confirmation process in this review.**

The code can check catalogs, fit a transit model, compare sectors and products,
search secondary eclipses, inspect Gaia neighbors and pixel centroids, compute
TRICERATOPS-style false-positive probabilities, predict future transits, and
produce a human-reviewed CTOI checklist. The gate explicitly labels human
steps as required; a passing automated report is still only a candidate.

The benchmark contains already-confirmed planets and no new candidate. It
contains no TPF or FFI pixels, no follow-up photometry/spectroscopy, and no
independent observer reproduction. The official [SPOC products](https://archive.stsci.edu/hlsp/tess-spoc)
include target pixels and, for detected TCEs, data-validation products. QLP
also publishes validation reports; a light curve from one product does not
replace those source-local tests.

TRICERATOPS is a Bayesian scenario model, not a generic confidence score.
Published use combines TESS aperture/neighbor information and can benefit from
high-resolution imaging or multicolor follow-up. No FPP value from the local
software alone establishes discovery.

## Promotion states

Use the repository's intended distinction:

1. **Detected signal:** automated statistic and candidate ephemeris.
2. **Vetted candidate:** catalog checks, independent sector/product, pixel and
   astrophysical false-positive checks documented.
3. **Independently validated candidate:** appropriate statistical validation
   plus an independent pipeline or observing evidence.
4. **Confirmed discovery:** human review and external evidence sufficient for
   the object class; submit and publish only with evidence provenance.

No current benchmark result is at stages 2–4.

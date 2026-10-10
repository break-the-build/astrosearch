# Discovery readiness

## Status

**Not ready to claim a discovery.** This review evaluated confirmed planets and
matched inverted controls; it did not identify and validate a new object.

## Current stage definitions

1. **Detected signal:** periodic dip/event plus search statistic and reproducible
   input record.
2. **Vetted candidate:** catalog checks; independent sector/pipeline agreement;
   odd/even, secondary, sector-depth and centroid/pixel checks; plausible stellar
   and blend analysis.
3. **Independently validated candidate:** statistical false-positive analysis
   with the needed aperture/neighborhood and follow-up constraints; independent
   reproduction or confirmation appropriate to the object class.
4. **Confirmed discovery:** external observational evidence, human disposition,
   documented provenance and publication/announcement process appropriate to
   the field.

The Rust candidate row is stage 1. `ctoi_gate.py` can assemble a stage-2/3
checklist, but its automated pass is not a human confirmation. An opt-in
SQLite registry now stores hashed runs, candidate measurements, embedded JSON
provenance, and append-only status history. It has only been exercised by
synthetic unit fixtures; no science candidate or benchmark object has been
imported. Use `scripts/candidate_registry.py ingest` with a per-target
provenance manifest before starting a multi-session campaign.

## Candidate record required

For every signal retain a stable candidate identifier, TIC/source coordinates
and frame, TESS sectors/products, input and output checksums, source revision,
search configuration, time convention, period/epoch/duration/depth and
uncertainties, statistic/tier, event evidence, alternative aliases, vetting
checks, catalog cross-matches with query date, pixel/neighbor evidence, FPP
assumptions, independent reproduction, disposition history, and follow-up
links. Never set “confirmed” solely from a classifier or score.

## Before a CTOI

Run the project gate, but verify its product and catalog status dates manually;
check a current TESS/SPOC/QLP path and pixel-level evidence; compute FPP only
with appropriate contaminant and imaging priors; publish the predicted transit
windows with a timestamp; then submit the CTOI and seek follow-up. A high score
or absence from one catalog is insufficient.

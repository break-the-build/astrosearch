# Astrosearch survey summary — 2026-10-10

## Scope and reproducibility

A 25-target QLP screen was selected from a TIC sample of nearby K/M dwarfs using actual MAST product coverage. Normal survey settings were used (Sector 90+ and multi-year baseline), followed by a separate SPOC 2-minute search for the seven QLP signal-bearing targets. The QLP batch completed 24 targets; one selected target failed the live recent-sector/baseline check and one eligible target was deferred by the 25-target limit. The run metadata, target ranking, coverage cache, per-target logs, and SQLite ledgers preserve provenance.

The release binary was `983934f86cb0a31d328bc9ab9c1382c013a6ec4d028f78b14ec1b570e5eeec5d`; the source tree was dirty during the run. Do not use this run as a version-to-version performance benchmark.

## What the searches recovered

The independent SPOC pass recovered known periodic signals on:

- TIC 150428135: TOI 700.01, about 16.051 d.
- TIC 101955023: GJ 1132 b, about 1.629 d.
- TIC 307210830: L98-59 b, c, and d, about 2.253, 3.691, and 7.451 d.
- TIC 413248763: GJ 357 b, about 3.931 d.

These are known-signal recoveries, not new objects. Four of the seven signal-bearing targets also have matching recent SPOC TCE catalog rows. The 1.049-day signal on L98-59 matches SPOC TCEs but is not in the confirmed-planet catalog, so it is not a new unrecorded signal.

## Unmatched signals requiring caution

| Target | Astrosearch signal | Independent/catalog checks | Assessment |
| --- | --- | --- | --- |
| TIC 272232401 | QLP 2.82418 d, rank 121; two adjacent candidate rows | SPOC has many nearby-period triggers (about 2.822–2.833 d); no matching recent SPOC TCE, TOI, or confirmed transiting planet | Not ready for follow-up as a planet candidate. The readiness gate finds highly inconsistent depths between sectors and QLP/SPOC, plus a strong odd/even difference. Its independent BLS check recovers a recurring dip, but that does not establish an astrophysical transit. Physical fit timed out; empirical FAP and pixel/FPP checks were not completed. |
| TIC 155315739 (GJ 887) | QLP 24.0109 d duo-transit, rank 25; also a dense 13–14 d cluster | No detection-tier SPOC recovery and no matching recent SPOC TCE | Weak follow-up lead only. The many 13–14 d rows are aliases/related triggers, not separate objects; several have extreme odd/even differences. |
| TIC 118585685 | QLP 19.787 d, only two events | No SPOC detection-tier recovery, TOI, or matching recent TCE | Reject as a planet lead for now: the odd/even difference is very large (74σ) and the signal has only two events. |

No signal in this survey passed the readiness checks for a credible new planet. The result is a useful recovery and vetting exercise, not a discovery claim.

## Code and operational fixes made

- Survey products now write `combined.meta.json` with per-sector and pipeline provenance, selected flux column, and crowding-correction status. Older combined files can infer sectors from their fifth CSV column.
- The candidate gate no longer counts the same QLP light curve twice when comparing QLP against SPOC.
- Transit-model fitting in the CTOI gate now has a configurable timeout and records timeout/nonzero exit as a structured fit failure, allowing later checks and a report to complete.
- Batch TCE preflight now retries the public index, reports match/no-match/partial coverage accurately, updates stale per-target notes, and preserves a prior complete catalog snapshot during an incomplete fetch.
- The MAST TCE parser recognizes comment-prefixed headers and the `tce_time0bt` epoch field.

The corrected TIC 272232401 gate is at `gates/TIC272232401-rank121-corrected/gate.md`. It was run offline for the automated gate portion, with no null trials or TRICERATOPS/pixel check; those are explicitly unevaluated. Separate public TOI/confirmed-planet and SPOC TCE preflight results are recorded in the run directory.

## Limits

The candidate threshold has limited known-planet recall in the existing benchmark (87/200 blind known planets recovered). A non-detection does not rule out a planet. The SPOC run used the updated 120-second SPOC threshold; the initial QLP screen predates that fix and was effectively QLP-only. Sector counts can be large and often include old sectors, so actual recent-sector and baseline eligibility must be checked per product stream. This transit search does not search generically for asteroids, supernovae, or other transient classes.

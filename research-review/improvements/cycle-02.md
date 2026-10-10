# Improvement cycle 02 — Persist candidate lineage and human dispositions

## Baseline and hypothesis

Search runs currently emit per-target candidate CSV files. The survey ledger
tracks target processing, but does not provide a stable cross-run candidate ID,
retain the observing manifest with the candidate, or record who advanced its
review status and on what evidence.

**Hypothesis:** an append-only SQLite registry with content hashes and gated,
human-attributed status transitions will preserve enough lineage for a
multi-session campaign without automatically turning detections into claims.

## Changes

- Added `scripts/candidate_registry.py` with `ingest`, `promote`, and `list`
  commands. It records input and candidate-output hashes, source revision,
  candidate measurements/evidence, and the JSON provenance manifest contents.
- Run and detection identifiers are deterministic for identical target,
  revision, input, candidate output, and candidate row. Re-ingesting the same
  run is idempotent.
- New rows begin at `DETECTED_SIGNAL`. Promotion advances exactly one stage
  and requires a human actor and an evidence reference. Benchmark rows were not
  imported and no candidate was promoted.
- Added synthetic unit coverage for idempotent ingest, manifest retention,
  required human evidence, and sequential stage transitions.

## Acceptance and result

- Three registry-specific unit tests pass. Full Python unit suite: 64 passed.
- This is an opt-in storage tool, not an automated validator or candidate
  classifier. It has not yet been integrated with the survey runner, and it
  has not been exercised against an actual search output. Adoption requires
  providing a provenance JSON manifest for each ingested run.
- The registry does not enforce scientific sufficiency of evidence text; a
  human can still enter an inadequate reference. Audit evidence and human
  dispositions in review.

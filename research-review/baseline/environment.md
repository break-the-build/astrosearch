# Environment and repository state

## Frozen source baseline

- Commit: `fc0e95e99faddb92f20d5b94d6319ce4df33c9fc`.
- Branch at inspection: `main`.
- Review date: 2026-10-10 (America/Los_Angeles).
- No tracked changes at inspection.
- Existing untracked `4206066_gaia_background.csv` preserved.
- `.env` and local generated artifacts are ignored; no credentials were read.

Creating a review branch failed because the workspace permission profile makes
`.git` read-only. No commit, push, merge, or branch switch was performed. New
review files and the small evaluator changes are in the writable working tree.

## Software versions

- macOS 26.6, arm64.
- Rust 1.98.1; Cargo 1.98.1.
- Python 3.14.7.
- Known-answer builder provenance records Lightkurve 2.6.0.
- The `.venv` contains the dependencies used for the Python unit tests and
  analysis; `pytest` itself is absent. The test suite uses `unittest` and ran
  without installing anything.

An attempted `uv run` dependency sync failed because DNS could not resolve
PyPI. The already-installed `.venv/bin/python` worked. Network-dependent archive
queries were not repeated.

## Hardware and capacity

`system_profiler` reports a Mac Studio with Apple M2 Ultra, 24 CPU cores
(16 performance, 8 efficiency), 60 GPU cores, and 128 GB RAM. The GPU has Metal
support, but the search code has no GPU implementation. At inspection, the
workspace volume had about 327 GB free; the generated benchmark occupies about
2.4 GB. No paid compute or external spend was used.

The benchmark search used 24 independent Rust processes with one Rust thread
per target, one transit pass, 0.15–30 day bounds, event grouping, and a
250,000-trial cap. Sidecars preserve binary, input, and per-process search
configuration hashes. They do not preserve per-target runtime or verify the
24-process concurrency setting; the latter is recorded in benchmark
provenance.

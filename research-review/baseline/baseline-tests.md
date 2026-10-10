# Baseline and post-change checks

## Before review changes

**VERIFIED.** `cargo test --release` passed 62 Rust tests: 39 library tests,
5 correlated-noise integration tests, 4 mixed-cadence tests, and 14 regression
tests. No failures or ignored tests.

The repository's Python tests are written with `unittest`. Before making review
changes, `.venv/bin/python -m unittest discover -s tests -p 'test_*.py' -v`
passed 57 tests. `pytest` is not installed, so `python -m pytest` failed with
`No module named pytest`; no package install was attempted after network/DNS
failure.

## After review changes

The evaluator changes and registry cycle passed 64 Python unit tests. In the
subsequent recommendation-implementation cycle, the modified and new scripts
passed `py_compile`; the full unit suite was not rerun. The Rust detector was
not modified by these Python-only changes, so the 62-test Rust baseline remains
applicable; `cargo test --release` was not rerun after Python-only changes.

The benchmark itself was rerun from the frozen search outputs with two
statistics. Evaluation scripts validate target-level ephemerides against the
truth manifest; full metrics are in the ignored release directory. The
benchmark output auditor checked all 600 required rows, all light-curve input
hashes, and all candidate-output hashes.

## What these checks establish

Passing unit and regression tests supports software behavior on the encoded
fixtures. It does not validate every TESS product, establish external
reproducibility, calibrate astrophysical false positives, or prove that the
search grid is sufficiently resolved for the full period-duration domain.

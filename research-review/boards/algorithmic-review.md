# Board 5 — Algorithm and signal-processing review

**Sequential simulated review; not an independent external audit.**

## Finding

**Rating: 2.5/5 — thoughtful signal-processing implementation, but the current
coarse grid under-resolves the declared search space for nearly all tested
targets.**

The implementation shares preprocessing across detectors, estimates noise by
segment and transit duration, groups ephemerides by shared observed events,
labels aliases, and can iteratively mask signals. The sinusoid path uses a
recurrence. These are sensible engineering choices and the Rust regression
suite checks serial/parallel determinism.

The observed grid cap is material: 583 of 600 required searches hit 250,000
trials; the median requested count was 4.28 million. The current period grid
tries to limit phase drift but caps the number of points; the resulting warning
explicitly says narrow peaks may be missed. The search is therefore not fully
resolving its own stated timing criterion. Better parallelism alone cannot fix
the missing period samples.

TLS is present as an optional independent path, not the Rust core. Published
TLS benchmarking reports better small-planet recovery under its test
conditions; that result does not transfer automatically to this TESS sample.
The local Wōtan-inspired biweight implementation likewise needs local
injection/recovery comparison.

## Recommendations

1. Replace the uniform 250,000-period sweep for long baselines with a validated
   coarse-to-fine or event-time-based search. Preserve a dense local refit for
   shortlisted peaks. Measure the full trial correction and alias behavior.
2. Use the cap as an explicit candidate-quality flag and report expected timing
   drift at each trial rather than silently treating capped output as complete.
3. Compare box and TLS on held-out targets at the same host-level null trigger
   rate. Include active stars, long baselines, gaps, and multi-planet cases.
4. Use the current known-answer set for descriptive recovery only; it is no
   longer blind after this review.

**HYPOTHESIS:** event-time or coarse-to-fine methods could search the full
declared period range faster while resolving narrow peaks. A paired benchmark
must establish both recall and false alarms before replacing the current grid.

Sources: [TLS paper](https://arxiv.org/abs/1901.02015),
[Wōtan benchmark](https://arxiv.org/abs/1906.00966).

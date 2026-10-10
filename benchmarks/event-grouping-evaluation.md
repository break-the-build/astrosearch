# Event-grouping paired evaluation

## Setup

Compared the existing ephemeris grouping mode with the new opt-in event grouping mode using the calibration harness's paired deterministic seeds:

```sh
cargo run --release --example correlated_noise_calibration -- --scale 0.02 --real-limit 4 --grouping ephemeris --out target/p0-3-final-ephemeris
cargo run --release --example correlated_noise_calibration -- --scale 0.02 --real-limit 4 --grouping events --out target/p0-3-final-events
```

Both modes saw the same white-noise, correlated-noise, two-event, two-planet, and real inverted-control cases. Each invocation reported 24 serial/parallel checks with zero mismatches. The calibration outputs are under the ignored `target/` directory.

## Results

- False-alarm null thresholds and counts were identical between grouping modes. The grouping option does not alter the score statistic.
- Recovery on white-noise injections was identical in the sampled cases.
- On correlated-noise multi-event injections, event grouping retained 8/12 at top 12 versus 8/12 for ephemeris grouping; at top 1 it was 5/12 versus 5/12. Long-duration correlated cases also tied: 5/12 recovered at top 12 for 4 h and 6 h injections in both modes.
- Real inverted-curve comparisons tied: 4/8 at 2 h and 1/8 at 4 h at the reported matched threshold.
- Total elapsed calibration time was 332 s for event grouping and 321 s for ephemeris grouping (+3.4%). The correlated-scoring searches took 99 s and 93 s respectively (+6.5%). These are single-run timings and are not a stable benchmark.

## Interpretation

This small run does not show an accuracy or false-alarm advantage for event grouping. It does show no regression in the aggregate recovery counts in the cases above, at a modest observed runtime cost. The sample is too small to claim equivalence or calibrate a production threshold: most strata contain 4–12 injections, so their Wilson intervals are wide.

Event grouping is therefore available as an explicit candidate-output mode, while the established ephemeris mode remains the default. Its main benefit is preserving event-level evidence and alternate ephemerides without consuming separate candidate slots; the synthetic suite does not establish improved planet yield. Before changing the default, expand injections across periods, transit phases, gaps, sector boundaries, signal depths, and same-period blended/multiple-event cases, then run on a larger set of real light curves.

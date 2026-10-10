# Partial search evaluation (2026-10-10)

Both long-running searches were interrupted at the user's request. Completed per-case outputs were retained and evaluated below. Neither report is a final benchmark result.

## TLS comparison

- Completed: 90 of 600 planned cases.
- All 90 completed cases are calibration inverted/null controls; no blind positive or blind control cases had completed.
- A provisional TLS SDE threshold of **45.081** triggers 4/90 controls (**4.4%**, Wilson 95% CI 1.7–10.9%).
- The frozen Astrosearch jackknife-score threshold triggers 4/90 of the same controls (**4.4%**).
- The paired controls split 3 TLS-only, 3 Astrosearch-only, 1 both, and 83 neither. The exact McNemar p-value is 1.0; the paired bootstrap 95% interval for the trigger-rate difference is −5.6 to +5.6 percentage points.
- This is a descriptive calibration-control check, not evidence of equal performance: TLS's threshold was selected on these same 90 controls, its SDE is not on the Astrosearch score scale, and there are no blind planet recoveries. Sensitivity and any benefit from TLS remain unmeasured.

Machine-readable details: [TLS partial evaluation JSON](generated-20261009/tls-full/evaluation-partial.json).

## Eclipsing-binary hard negatives

- Completed: 12 of 136 targets with supported light-curve products. The broader catalog selection had 200 targets; 64 lacked a supported product and were excluded before searching.
- At the frozen Astrosearch threshold, **12/12** known eclipsing binaries triggered the score threshold (**100%**, Wilson 95% CI 75.8–100%).
- This is expected for a periodic-dip search: it shows that the search stage finds eclipses. It also means the score threshold alone cannot distinguish planets from eclipsing binaries. The next vetting stages must reject them using secondary eclipses, odd/even depths, centroid/blend checks, and other evidence.
- The sample is small and incomplete, so this rate is preliminary and should not be treated as catalog-wide performance or as a post-vetting false-positive rate.

Machine-readable details: [hard-negative partial evaluation JSON](hard-negatives/partial-20261010/evaluation-partial.json). Frozen labels and hashes are in [the partial truth manifest](hard-negatives/partial-20261010/truth-manifest.csv).

## Run state

Both parent jobs were stopped with Ctrl-C; their completed outputs are resumable. The TLS output directory contains per-case raw results and checksummed sidecars. The hard-negative survey contains per-target manifests and independent chunk ledgers. Resume with the same scripts and output paths to reuse completed work.

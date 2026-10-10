# Recommendation implementation status

| Review recommendation | Status in this working tree | Evidence still needed |
| --- | --- | --- |
| Period-grid cap | Quantified; remains capped and reported in warnings | Fresh injection and dense-oracle comparison before accepting a new algorithm |
| Hard astrophysical negatives | Manifest freezer and per-class metrics implemented | Actual externally dispositioned EB/blend/variable/artifact sample |
| Independent search methods | TLS known-answer runner added; QLP/SPOC survey and TLS paths are available | Run TLS on a frozen common set; source both pipeline products for matched hosts |
| Exact release identity | New-run source fingerprint, binary hash, compiler, platform and package versions recorded | Historical binary/source association and external immutable archive/DOI |
| Per-target metrics | New Rust searches record wall time, sampled CPU/RSS, config and exit status | Re-run on new inputs; old release has no retroactive timing data |
| QLP product-era strata | Builder tags product family; evaluator can backfill from case metadata | A pure S94+ host sample and matched re-extraction are absent |
| Injection/recovery | Existing real-curve harness samples period, phase, duration, depth, impact parameter, exposure, gaps and season holdout | Preregistered multi-target grid and uncertainty summaries |
| Candidate ranking/depth uncertainty | Evaluator reports top-1/5/10/any-rank; depth uncertainty still not calibrated | Fresh holdout with physical-depth truth and sector-level reference estimates |
| Physical vetting | Synthetic checks exist | Blinded archival dispositions and reviewer agreement |
| GPU/Metal | Deferred | Verified search-resolution redesign and CPU-oracle equivalence |

## Descriptive rerun completed

`evaluate_benchmark.py` was rerun on the frozen outputs with
`case-metadata.csv`. It preserves the calibration-only threshold and reports
candidate rank plus marginal product-era, cadence, baseline, sector-count,
period, depth, and magnitude strata. Results are in the ignored local output
`data/benchmarks/generated-20261009/evaluation-jackknife-by-product.json` and
are summarized in `../experiments/benchmark-results.md`.

This is a re-analysis of the already inspected release. It does not restore a
blind holdout or validate any new threshold.

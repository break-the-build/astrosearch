# CTOI gate: TIC 272232401

Candidate: periodic_dip rank 121, P = 2.824181 d, score 46.10

**Verdict:** Not ready for a CTOI: 1_multisector_multiproduct, 2_physical_fit_density, 3_no_false_positive_indicators, 4_not_already_known, 5_significance, 6_preregistered_ephemeris.

> Automated screening only. Passing every check does not make this a planet; it makes it worth a human vetting review and follow-up.

| # | Criterion | Status |
| --- | --- | --- |
| 1 | Signal in 2+ sectors and 2+ products, consistent depth | FAIL |
| 2 | Physical transit fit preferred; stellar density consistent | NOT_EVALUATED |
| 3 | No secondary, odd/even, grazing, centroid, neighbor, RUWE flags | FAIL |
| 4 | Not a known TOI, CTOI, TCE, or confirmed planet | NOT_EVALUATED |
| 5 | Empirical FAP and TRICERATOPS FPP below 1% | NOT_EVALUATED |
| 6 | Ephemeris with uncertainties, predicted windows | FAIL |

## Sub-checks

- 1_multisector_multiproduct / two_sectors_with_events: **PASS**
- 1_multisector_multiproduct / depth_consistent_across_sectors: **FAIL**
- 1_multisector_multiproduct / two_independent_products: **FAIL**
- 2_physical_fit_density / transit_model_preferred_delta_bic: **NOT_EVALUATED**
- 2_physical_fit_density / stellar_density_consistent: **NOT_EVALUATED**
- 3_no_false_positive_indicators / no_secondary_eclipse: **PASS**
- 3_no_false_positive_indicators / odd_even_consistent: **FAIL**
- 3_no_false_positive_indicators / not_grazing: **NOT_EVALUATED**
- 3_no_false_positive_indicators / not_dominated_by_momentum_dumps_or_scattered_light: **PASS**
- 3_no_false_positive_indicators / no_centroid_offset: **NOT_EVALUATED**
- 3_no_false_positive_indicators / no_neighbor_able_to_mimic_depth: **PASS**
- 3_no_false_positive_indicators / ruwe_below_1_4: **NOT_EVALUATED**
- 5_significance / empirical_fap_upper95_below_1pct: **NOT_EVALUATED**
- 5_significance / triceratops_fpp_below_1pct: **NOT_EVALUATED**
- 5_significance / score_at_or_above_calibrated_detection_tier: **PASS**

## Derived planet parameters (if on target)


## Verification protocol (docs/VERIFICATION_PROTOCOL.md)

| Step | Requirement | Status |
| --- | --- | --- |
| 1 | Provenance: committed code revision and input SHA-256 recorded | FAIL |
| 2 | Score at or above the calibrated detection tier | PASS |
| 3 | All automated gate criteria pass | FAIL |
| 4 | Independent algorithm (BLS) recovers the signal | PASS |
| 5 | Difference-image centroid on target | NOT_EVALUATED |
| 6 | Independent reproduction from MAST by a second person | ACTION_REQUIRED |
| 7 | Timestamped public ephemeris (e.g. Zenodo) before follow-up | ACTION_REQUIRED |
| 8 | Human ExoFOP CTOI submission and TFOP follow-up | ACTION_REQUIRED |

## CTOI fields (draft, for human review)

- tic_id: 272232401
- period_days: 2.82418051
- period_unc_days: None
- epoch_bjd: 2458412.23449315
- epoch_unc_days: None
- depth_ppm: 9398.464810991629
- depth_unc_ppm: 5.195270447379306
- duration_hours: 8.14158288
- planet_radius_earth: None
- planet_radius_unc_earth: None
- notes: astrosearch periodic_dip

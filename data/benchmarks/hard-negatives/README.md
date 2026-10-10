# TESS eclipsing-binary hard negatives

This cohort is drawn from table 3 of Kostov et al. (2025), the TESS Ten
Thousand Catalog: 7,936 new eclipsing binaries that passed the paper's
ephemeris, photocenter, and visual-vetting steps. The VizieR response is frozen
as `tess-eb-catalog-table3.tsv`; its SHA-256 and selection rules are recorded
in `selection.json`.

The build is reproducible in stages:

1. Download the public table from VizieR:

   ```sh
   curl -L --fail 'https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=J/ApJS/279/50/table3&-out.max=10000' \
     -o data/benchmarks/hard-negatives/tess-eb-catalog-table3.tsv
   ```

2. Make a deterministic sample (the current design selects 200 targets over
   period and TESS-magnitude strata), excluding every host in the frozen TESS
   planet catalog and the known-answer sample:

   ```sh
   .venv/bin/python scripts/prepare_tess_eb_hard_negatives.py --limit 200
   ```

3. Query MAST once per TIC to select an available QLP, TESS-SPOC, or SPOC
   reduction, preferring QLP and recording targets without supported products:

   ```sh
   .venv/bin/python scripts/probe_tess_lightcurve_products.py \
     --targets data/benchmarks/hard-negatives/survey-targets.csv \
     --output data/benchmarks/hard-negatives/probed-eb-targets.csv --workers 6
   ```

4. Run the selected reductions through the same Rust search used by the survey:

   ```sh
   .venv/bin/python scripts/run_hard_negative_survey.py \
     --targets data/benchmarks/hard-negatives/probed-eb-targets.csv \
     --output data/benchmarks/hard-negatives/survey-parallel --workers 4 --threads 6
   ```

5. Freeze successful outputs and evaluate with a threshold calibrated on
   separate known-answer inverted controls using the same binary and search
   settings:

   ```sh
   .venv/bin/python scripts/build_survey_hard_negative_labels.py \
     --survey-output data/benchmarks/hard-negatives/survey-parallel \
     --selection data/benchmarks/hard-negatives/probed-eb-targets.csv \
     --output data/benchmarks/hard-negatives/labels.csv
   .venv/bin/python scripts/build_hard_negative_manifest.py \
     --labels data/benchmarks/hard-negatives/labels.csv \
     --output data/benchmarks/hard-negatives/truth-manifest.csv
   .venv/bin/python scripts/evaluate_benchmark.py \
     --manifest data/benchmarks/hard-negatives/truth-manifest.csv \
     --fixed-threshold THRESHOLD_FROM_MATCHING_CALIBRATION \
     --output data/benchmarks/hard-negatives/evaluation.json
   ```

The hard-negative trigger rate measures how often an eclipsing binary produces
a candidate above the planet-search screen. It is not the false-positive rate
of a complete candidate vetter; the detector is expected to find eclipses, and
secondary-eclipse, odd/even, centroid, and catalog checks must reject them.
Keep `survey.log`, target manifests, hashes, and per-case outputs with the
evaluation. Do not combine products from different pipelines into one light
curve; the selected author is recorded in each target manifest.

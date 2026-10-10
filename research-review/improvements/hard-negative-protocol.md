# Hard-negative cohort protocol

`evaluate_benchmark.py` already separates `known_false_positive` outcomes from
null triggers and now reports trigger rates by `false_positive_type`. Use
`build_hard_negative_manifest.py` to freeze externally curated labels only
after the corresponding curve and astrosearch output exist.

The labels CSV must contain:

```text
case_id,host_id,split,object_type,label_reference,data_csv,candidates_csv
```

Use `split=calibration` or `split=blind` at the host level. `object_type`
should identify a documented class such as eclipsing binary, blend, variable,
or instrumental artifact. `label_reference` must point to the catalog or
archived disposition that supports that label. The script records input and
candidate-output SHA-256 hashes and rejects reused hosts across splits. It
does not infer or adjudicate a label.

Example:

```sh
python scripts/build_hard_negative_manifest.py \
  --labels data/benchmarks/hard-negatives/labels.csv \
  --output data/benchmarks/hard-negatives/truth-manifest.csv

python scripts/evaluate_benchmark.py \
  --manifest data/benchmarks/hard-negatives/truth-manifest.csv \
  --fixed-threshold 22.48764564 \
  --output data/benchmarks/hard-negatives/evaluation.json
```

For that standalone command, the fixed threshold must have been selected on a
separate calibration-null cohort with the same statistic and search settings.
Alternatively, combine the false-positive rows with a compatible known-planet
and null manifest before evaluating. The current local benchmark contains no
externally dispositioned false-positive cases, so no real hard-negative rate
is reported yet.

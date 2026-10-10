# Board 2 — Software engineering review

**Sequential simulated review; not an independent external audit.**

## Finding

**Rating: 3.5/5 — solid prototype engineering with an unresolved scale/accuracy
constraint.**

**VERIFIED:** Rust code separates ingestion, preprocessing, segmentation,
period-grid construction, transit search, variability, grouping, refinement,
provenance, and output. It has safe path-collision checks, finite-value
filtering, deterministic serial/parallel regression checks, resumable Python
survey operations, atomic output handling, and 62 passing Rust tests. The
survey runner stores sector products, manifests and checksums. The 600-case
benchmark is resumable and hashes inputs, executable, configuration, and
outputs.

**VERIFIED DEFECT:** 97.2% of the benchmark searches requested 599,000–6,558,648
period trials but were capped at 250,000. The sidecar does not record per-target
runtime; the benchmark manifest also does not associate the binary hash with a
source commit or preserve the executable. This weakens exact rerun
reproducibility.

**VERIFIED OPERATIONAL LIMIT:** the repository's Git directory is read-only in
this execution, so the branch requested by the local review protocol could not
be created. The starting commit and untracked user file were recorded and
preserved.

## Recommendations

- Have the runner store elapsed wall/CPU time, peak RSS if available, exit
  status, warning counts, and machine/runtime version in each sidecar.
- Record a source commit and build flags alongside the binary hash; archive the
  exact executable for immutable benchmark releases.
- Keep one output schema version and make the evaluator reject missing or
  malformed per-case outputs rather than silently shrinking denominators.
- Add a benchmark report that checks cap frequency and summarizes throughput
  without merging across substantially different cadences/baselines.
- Do not increase worker count beyond measured CPU/memory limits. This workload
  is CPU-bound; the Metal GPU is unused.

No correctness regression was observed in the existing suite. The lack of
runtime telemetry is an observability gap, not proof that execution is
unreliable.

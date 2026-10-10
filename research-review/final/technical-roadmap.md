# Technical roadmap

## Next release sequence

### 1. Fix search-space resolution first

Design a full-resolution coarse-to-fine or event-based period search. Define the
period-duration domain and timing error bound explicitly. Keep a dense local
refinement/oracle and emit a warning/quality flag whenever that bound is not
met. Use fresh injection phases and targets; do not tune against the retired
known-answer holdout.

### 2. Build an astrophysical hard-negative set

Add catalogued eclipsing binaries, known blends, variable stars and archived
false alarms. Keep signal detection and planet classification as separate
metrics. Compute per-class rates and host-clustered confidence intervals.

### 3. Compare independent search paths

Run Rust box search, optional TLS, and archived QLP/SPOC TCE products on the
same normalized inputs. Calibrate each method on separate controls. Compare
paired recovery and false triggers; document search parameters and pipeline
versions.

### 4. Improve data/operational provenance

Store source commit, executable, compiler flags, CPU thread/worker settings,
wall time, CPU time, peak memory, trial demands, warnings, sector/product
versions, time system, flux column and quality mask in each run record. Archive
an immutable release and copy its checksums to the report.

### 5. Harden and adopt the candidate lifecycle registry

Create stable candidate IDs and append-only evidence/status history across
surveys. Distinguish detected signal, vetted candidate, independently
validated candidate and confirmed discovery. Keep manual disposition fields
human-owned; candidate rows are not automatically promoted.

### 6. Recalibrate and validate

After algorithm changes, run injections across period, duration, depth,
impact parameter, cadence, phase, gaps and stellar variability. Use a new blind
known-host split. Report completeness by target class and use confidence
intervals. Freeze this benchmark release as baseline v1.

### 7. Follow up only evidence-backed candidates

Use catalog and pixel checks, secondary/odd-even tests, sector consistency,
neighbor analysis, appropriate FPP constraints, independent search, timestamped
predictions, and human observing review. Preserve negative results.

GPU/Metal evaluation comes after period-grid redesign. The hardware GPU is
available, but current output accuracy is limited by missing trial samples;
accelerating the capped search does not restore those samples.

## Improvement-cycle record

Three implementation cycles are documented in
[cycle-01.md](../improvements/cycle-01.md),
[cycle-02.md](../improvements/cycle-02.md), and
[cycle-03.md](../improvements/cycle-03.md). The third cycle adds workflow and
evaluation tools but has not produced a new TLS cohort result or a real
hard-negative sample. Further scientific algorithm cycles should use a new
hypothesis, preregistered acceptance criteria, and fresh blind or injection
data. The ten-cycle target is not declared complete.

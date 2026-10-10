# Verification protocol for an astrosearch candidate

astrosearch finds **candidates**. None of its outputs is a planet. This
protocol lists what has to be true before anyone describes a signal as a
community candidate (an ExoFOP CTOI). It also lists the human steps that come
after that. `scripts/ctoi_gate.py` fills in steps 1–5 automatically: it writes
them to `gate.json` (`verification_protocol`) and to `gate.md`. Steps 6–8
always stay `ACTION_REQUIRED`.

| Step | Requirement | How it is checked | Pass condition |
| --- | --- | --- | --- |
| 1 | Provenance | The search writes `<candidates>.provenance.json` (git commit, input/config SHA-256). The gate records the commit, input digests and package versions. | Committed, clean tree (`git_dirty: false`) |
| 2 | Calibrated significance | `tier` column in the search CSV; gate criterion 5 | Jackknife score ≥ detection tier (`data/calibration/detection_threshold.json`, currently 23.3) |
| 3 | Automated vetting | The six gate criteria | All automated criteria `PASS` |
| 4 | Independent algorithm | astropy BLS cross-check in the gate | Candidate period (or P/2, 2P) among the top 5 BLS peaks |
| 5 | On-target signal | Difference-image centroid on up to 2 SPOC sectors (`--pixel-sectors`) | No offset that is both larger than 1 pixel and ≥ 3σ |
| 6 | Independent reproduction | A second person re-downloads from MAST, checks the SHA-256 digests against the provenance records and re-runs the recorded commands | Same candidate, same tier, same gate verdict |
| 7 | Pre-registration | Publish the ephemeris and predicted windows (`ephemeris_windows.csv`) with a timestamp, for example on Zenodo or OSF | Public, timestamped record before any follow-up |
| 8 | Community follow-up | Submit the CTOI on ExoFOP yourself (never automated), then follow up through TFOP: seeing-limited photometry, reconnaissance spectroscopy, high-resolution imaging | TFOP disposition |

## Commands

```bash
# 1–2: search with provenance; the tier column marks detection vs screening
cargo build --release
target/release/astrosearch --input data/<target>/combined.csv --output out/candidates.csv

# 3–5: gate with MCMC uncertainties, 300 nulls, FPP and pixel checks (network)
python scripts/ctoi_gate.py --tic <TIC> --lightcurve data/<target>/combined.csv \
    --candidates out/candidates.csv --rank <rank> --out-dir out/gate

# 6: reproduction: compare digests, then repeat both commands at the recorded commit
git checkout <git_commit from gate.json>
shasum -a 256 data/<target>/combined.csv   # must equal provenance.inputs.lightcurve.sha256
```

## What the checks do and do not establish

- **Detection tier.** The statistic is the jackknife score: the score with
  the strongest single event removed. It was calibrated on the known-answer
  benchmark (QLP light curves, one search pass). With the threshold fixed on
  the calibration split, 7.0% of blind-split null controls (14/200) and 43%
  of blind confirmed planets (86/200) reached it. The deployed threshold
  (23.3) gives a 5% trigger rate across all 400 controls. The tier controls
  false alarms from systematics. It does not say whether a real dip comes
  from a planet.
- **Secondary-eclipse scan.** It wrongly fails 2.5% of blind benchmark
  planets (5/199) and catches 83% of injected half-depth secondaries
  (`scripts/evaluate_secondary_scan.py`). It cannot see a secondary that is
  shallower than the noise.
- **Blend tests.** These cover catalogued TIC neighbors only, using the flux
  ratio and the Vanderburg et al. (2019) shape limit. They cannot rule out
  uncatalogued companions, which is why steps 5 and 8 exist.
- **TRICERATOPS FPP.** This is a statistical estimate under its own priors.
  It does not validate a planet.

Passing steps 1–5 makes a signal worth a human vetting review, and nothing
more.

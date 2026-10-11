# Astrosearch agent instructions

- Start with `docs/GETTING_STARTED.md` and `docs/VERIFICATION_PROTOCOL.md`; consult `README.md` for implementation details and `docs/SURVEY_ORCHESTRATION.md` for batch/ledger behavior.
- Astrosearch finds and ranks light-curve signals. Call unconfirmed results **candidates/signals**, never planets. Never submit to ExoFOP. Never commit `.env`, credentials, or tokens.
- The recommended default is Profile J. Preserve the exact product, binary, settings, pass count, and calibration contract. A tier is calibrated only when those match; do not silently weaken contract checks or tune against the blind split.
- Keep data provenance, measured usable rows, null/negative controls, injection completeness, source localization, and independent predictive follow-up explicit. Metadata coverage is not usable data. Record missing evidence, errors, and timeouts.
- Survey runs use a persistent ledger and cache. Inspect them before starting work; reuse the same ledger/output to resume and avoid duplicate searches. Do not mark incomplete work complete.
- For algorithm changes, evaluate the end-to-end pipeline on target-grouped splits and realistic positives plus eclipsing-binary, blend, and systematic negatives. Compare at matched false-alarm rates and report confidence intervals, paired gains/losses, and measured compute cost. Keep calibration artifacts versioned and write new results separately until reviewed.
- Do not run long surveys, overwrite benchmark/calibration outputs, or change thresholds based on the blind split without explicit user direction. Do not fabricate counts, candidate parameters, runtime, or significance.
- Prefer small pilots, preserved logs/manifests, resumable bounded work, and concise documentation updates alongside user-facing behavior changes.

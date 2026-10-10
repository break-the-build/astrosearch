//! Candidate screening for positive-flux astronomical light curves.
//!
//! Pipeline: `input` reads and cleans samples; `preprocess` segments them
//! (`segments`), detrends, weights them, and estimates correlated noise
//! (`rednoise`). Detectors then run on the prepared samples: `transit`
//! (box search over the `grid` of trial periods, with `epochs` vetting and
//! `grouping` deduplication), `variability` (sinusoids), and `single_events`.
//! `refine` remeasures transit depths, and `output` renders CSV.

mod config;
mod epochs;
mod grid;
mod grouping;
mod input;
mod math;
mod output;
mod parallel;
mod phase;
mod preprocess;
mod provenance;
mod rednoise;
mod refine;
mod segments;
mod single_events;
mod transit;
mod types;
mod variability;

pub use config::{
    CandidateGrouping, Config, Scoring, TrendMethod, DETECTION_THRESHOLD, MAX_PERIOD_TRIALS,
};
pub use output::{render_csv, render_csv_with_tier};
pub use provenance::{provenance_json, sha256_hex};
pub use types::{
    AliasEphemeris, AliasType, Candidate, CandidateKind, EventEvidence, Point, Report,
    SegmentSummary,
};

/// Read a CSV or TSV light curve and perform candidate screening.
pub fn analyze_file(path: &str, config: Config) -> Result<Report, String> {
    config.validate()?;
    let (points, skipped) = input::read_points(path)?;
    let mut report = analyze(points, config)?;
    if skipped > 0 {
        report.warnings.push(format!(
            "skipped {skipped} rows with invalid values, nonpositive flux, or nonzero/invalid quality"
        ));
    }
    Ok(report)
}

/// Analyze one light curve. Repeated timestamps are conservatively median-combined.
pub fn analyze(points: Vec<Point>, config: Config) -> Result<Report, String> {
    config.validate()?;
    let original_count = points.len();
    let mut points = input::clean_points(points);
    if points.len() < 30 {
        return Err(format!(
            "need at least 30 valid, distinct observation times; found {}",
            points.len()
        ));
    }
    if config.invert {
        let center = math::median(points.iter().map(|point| point.flux));
        for point in &mut points {
            point.flux = (2.0 * center - point.flux).max(center * 1e-9);
        }
    }
    let span = points.last().unwrap().time - points[0].time;
    if !span.is_finite() || span <= 0.0 {
        return Err("observation times must have a finite, positive span".into());
    }
    let samples = preprocess::prepare_with_gap(
        &points,
        config.trend_days,
        config.segment_gap_days,
        config.trend_method,
    )?;
    if span / samples.cadence >= (1_u64 << 53) as f64 {
        return Err("time range is too large relative to the observation cadence".into());
    }
    let mut warnings = Vec::new();
    if points.len() != original_count {
        warnings.push(format!(
            "removed {} invalid rows or duplicate observations",
            original_count - points.len()
        ));
    }
    let transits = transit::iterative_periodic_dips(&samples, config, &mut warnings);
    let mut candidates = transits.candidates;
    refine::remeasure_depths(&mut candidates, &points, &samples, config);
    candidates.extend(variability::periodic_variability(
        &samples,
        config,
        &mut warnings,
    ));
    candidates.extend(single_events::events(&samples, config.min_event_points));
    single_events::suppress_periodic_event_duplicates(&mut candidates, samples.cadence);
    candidates.retain(|c| {
        [
            c.score,
            c.white_score,
            c.period,
            c.epoch,
            c.duration,
            c.depth,
            c.amplitude,
        ]
        .iter()
        .all(|v| v.is_finite())
            && c.odd_even_sigma.is_none_or(|v| v.is_finite())
            && c.noise_factor.is_none_or(|v| v.is_finite())
    });
    candidates.sort_by(|a, b| b.score.total_cmp(&a.score));
    let red = &samples.red_noise;
    let segments = samples
        .segments
        .iter()
        .zip(&samples.transit_noise)
        .enumerate()
        .map(|(index, (s, &noise_ppm))| SegmentSummary {
            start_time: points[s.start].time,
            end_time: points[s.end - 1].time,
            observations: s.len(),
            cadence_days: s.cadence,
            noise_ppm,
            red_noise: red.knots.iter().map(|&d| (d, red.beta(index, d))).collect(),
        })
        .collect();
    Ok(Report {
        observations: points.len(),
        span_days: span,
        transit_trials: transits.trials,
        segments,
        candidates,
        periodogram: transits.periodogram,
        warnings,
    })
}

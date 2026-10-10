//! Transit depths remeasured after re-detrending with each candidate's events masked.

use crate::{math, preprocess, preprocess::Samples, Candidate, Config, Point};

/// The box-fit depth comes from a trend fitted through the transits, which
/// can partly absorb them. For each periodic dip, re-detrend with its event
/// windows excluded and take the median of nearby out-of-transit samples (1.5-3
/// durations from an event) minus the median in-transit value. The result
/// replaces the box depth only with at least 2 in-transit and 4 nearby samples
/// and a positive value; those candidates say so in `note`.
pub(crate) fn remeasure_depths(
    candidates: &mut [Candidate],
    points: &[Point],
    samples: &Samples,
    config: Config,
) {
    for candidate in candidates.iter_mut().filter(|c| c.kind.is_transit()) {
        let excluded: Vec<bool> = points
            .iter()
            .map(|point| {
                candidate.event_evidence.iter().any(|event| {
                    (point.time - event.epoch).abs() <= candidate.duration / 2.0 + samples.cadence
                })
            })
            .collect();
        if !excluded.iter().any(|&value| value) {
            continue;
        }
        let detrended = preprocess::detrend_masked(
            points,
            samples,
            config.trend_days,
            config.trend_method,
            &excluded,
        );
        let mut inside = Vec::new();
        let mut outside = Vec::new();
        for (index, &value) in detrended.iter().enumerate() {
            let distance = candidate
                .event_evidence
                .iter()
                .map(|event| (points[index].time - event.epoch).abs())
                .fold(f64::INFINITY, f64::min);
            if distance <= candidate.duration / 2.0 {
                inside.push(value);
            } else if (1.5 * candidate.duration..=3.0 * candidate.duration).contains(&distance) {
                outside.push(value);
            }
        }
        if inside.len() >= 2 && outside.len() >= 4 {
            let depth = math::median(outside) - math::median(inside);
            if depth.is_finite() && depth > 0.0 {
                candidate.depth = depth;
                candidate.amplitude = depth;
                // Rows without this note keep the box-fit depth behind `score`.
                candidate
                    .note
                    .push_str("; depth remeasured after transit-masked detrending");
            }
        }
    }
}

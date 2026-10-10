//! Weighted floating-mean sinusoid search for periodic variability.

use crate::{
    grid::trial_count, grouping::collapse_periodic, math::robust_sigma,
    parallel::collect_trial_blocks, preprocess::Samples, Candidate, CandidateKind, Config,
};

pub(crate) fn periodic_variability(
    samples: &Samples,
    config: Config,
    warnings: &mut Vec<String>,
) -> Vec<Candidate> {
    let fmin = 1.0 / config.max_period;
    let fmax = (1.0 / config.min_period).min(0.5 / samples.cadence);
    if fmax <= fmin {
        return Vec::new();
    }
    // Five frequency samples per 1/baseline resolution element.
    let required = (5.0 * samples.span * (fmax - fmin)).ceil() as usize;
    let trials = trial_count(
        required.saturating_add(1),
        crate::grid::TrialBudget::from_config(&config),
        "variability",
        warnings,
    );
    let mut centered = samples.variability.clone();
    // Sector/pipeline offsets produce broad Fourier peaks that look periodic
    // across the stitched baseline. Remove each segment's weighted constant
    // while retaining variability inside each observing segment.
    for segment in &samples.segments {
        let range = segment.start..segment.end;
        let (mut sum, mut weight) = (0.0, 0.0);
        for (value, w) in centered[range.clone()]
            .iter()
            .zip(&samples.variability_weights[range])
        {
            sum += value * w;
            weight += w;
        }
        if weight > 0.0 {
            let offset = sum / weight;
            for value in &mut centered[segment.start..segment.end] {
                *value -= offset;
            }
        }
    }
    let noise = robust_sigma(&centered).max(1e-12);
    let weights: Vec<f64> = samples
        .variability
        .iter()
        .enumerate()
        .zip(&samples.variability_weights)
        .map(|((i, _), w)| {
            let z = centered[i].abs() / noise;
            w * if z > 4.0 { 4.0 / z } else { 1.0 }
        })
        .collect();
    let total_weight: f64 = weights.iter().sum();
    if total_weight <= 0.0 || !total_weight.is_finite() {
        return Vec::new();
    }
    let weighted_y: Vec<f64> = centered.iter().zip(&weights).map(|(y, w)| y * w).collect();
    let frequency_step = (fmax - fmin) / (trials - 1) as f64;
    let (step_sin, step_cos): (Vec<f64>, Vec<f64>) = samples
        .times
        .iter()
        .map(|time| (std::f64::consts::TAU * (frequency_step * time).rem_euclid(1.0)).sin_cos())
        .unzip();
    let candidates = collect_trial_blocks(
        trials,
        64,
        config.threads,
        || {
            (
                vec![0.0; samples.times.len()],
                vec![0.0; samples.times.len()],
            )
        },
        |(sines, cosines), block| {
            let first_frequency = fmin + frequency_step * block.start as f64;
            for (i, time) in samples.times.iter().enumerate() {
                let angle = std::f64::consts::TAU * (first_frequency * time).rem_euclid(1.0);
                (sines[i], cosines[i]) = angle.sin_cos();
            }
            let mut output = Vec::new();
            for k in block.clone() {
                let frequency = fmin + frequency_step * k as f64;
                let (mut cc, mut ss, mut cs, mut yc, mut ys, mut csum, mut ssum) =
                    (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
                for i in 0..samples.times.len() {
                    let (sin, cos) = (sines[i], cosines[i]);
                    let weight = weights[i];
                    cc += weight * cos * cos;
                    ss += weight * sin * sin;
                    cs += weight * cos * sin;
                    yc += weighted_y[i] * cos;
                    ys += weighted_y[i] * sin;
                    csum += weight * cos;
                    ssum += weight * sin;
                }
                let cc0 = cc - csum * csum / total_weight;
                let ss0 = ss - ssum * ssum / total_weight;
                let cs0 = cs - csum * ssum / total_weight;
                let determinant = cc0 * ss0 - cs0 * cs0;
                // Reject nearly singular projections, not just negative determinants.
                let candidate = if cc0 <= 0.0
                    || ss0 <= 0.0
                    || !determinant.is_finite()
                    || determinant <= 1e-12 * cc0 * ss0
                {
                    None
                } else {
                    let a = (yc * ss0 - ys * cs0) / determinant;
                    let b = (ys * cc0 - yc * cs0) / determinant;
                    let amplitude = a.hypot(b) * 1e6;
                    let score = (a * yc + b * ys).max(0.0).sqrt();
                    if score > 6.0 && score.is_finite() && amplitude.is_finite() {
                        Some(Candidate {
                    kind: CandidateKind::PeriodicVariability,
                    score,
                    white_score: score,
                    noise_factor: None,
                    period: 1.0 / frequency,
                    epoch: samples.origin,
                    duration: 0.0,
                    depth: 0.0,
                    amplitude,
                    n_events: (samples.span * frequency).floor() as usize,
                    event_evidence: Vec::new(),
                    alias_ephemerides: Vec::new(),
                    odd_even_sigma: None,
                    note: "sinusoidal screen after per-segment offset removal; check aliases and Gaia variability".into(),
                    })
                    } else {
                        None
                    }
                };
                if let Some(candidate) = candidate {
                    output.push((k, candidate));
                }
                for i in 0..samples.times.len() {
                    let sin = sines[i];
                    let cos = cosines[i];
                    sines[i] = sin * step_cos[i] + cos * step_sin[i];
                    cosines[i] = cos * step_cos[i] - sin * step_sin[i];
                    if (k - block.start + 1).is_multiple_of(32) {
                        let norm = sines[i].hypot(cosines[i]);
                        if norm > 0.0 {
                            sines[i] /= norm;
                            cosines[i] /= norm;
                        }
                    }
                }
            }
            output
        },
    );
    collapse_periodic(candidates, samples.span, 20)
}

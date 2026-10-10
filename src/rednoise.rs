//! Duration-dependent correlated ("red") noise for each segment.
//!
//! For bins as long as a transit, white noise predicts that a weighted bin mean
//! has variance 1/sum(w). Correlated noise inflates that by beta^2 (Pont, Zucker
//! & Queloz 2006). Beta is measured per segment at log-spaced duration knots and
//! interpolated linearly in log duration, so the transit search can apply a
//! segment- and duration-specific variance to every phase window.

use crate::{
    math::{median, robust_sigma},
    segments::{ratio, Segment, CADENCE_RATIO},
};

/// Longest searched transit duration, in days.
pub(crate) const MAX_TRANSIT_DURATION: f64 = 0.35;
const KNOTS: usize = 8;
/// Bin phase offsets; pooling them steadies the estimate for long durations.
const OFFSETS: usize = 4;
/// Fewer non-overlapping bins than this borrow a pooled estimate.
const MIN_BINS: usize = 10;
const MIN_BIN_POINTS: usize = 2;
/// Beta is never below 1 (no reward) and is capped at 10.
const MAX_BETA2: f64 = 100.0;

#[derive(Clone)]
pub(crate) struct RedNoise {
    /// Duration knots in days, increasing; one knot when only one duration exists.
    pub knots: Vec<f64>,
    /// beta^2, segment-major: `beta2[segment * knots.len() + knot]`.
    pub beta2: Vec<f64>,
}

impl RedNoise {
    pub fn estimate(
        times: &[f64],
        values: &[f64],
        weights: &[f64],
        segments: &[Segment],
        min_duration: f64,
    ) -> Self {
        let knots = knots(min_duration);
        // Normalized bin deviations for every (segment, knot), and how many
        // independent bins produced them.
        let samples: Vec<Vec<(Vec<f64>, usize)>> = segments
            .iter()
            .map(|s| {
                knots
                    .iter()
                    .map(|&width| deviations(times, values, weights, s, width))
                    .collect()
            })
            .collect();
        let pooled = |accept: &dyn Fn(&Segment) -> bool, knot: usize| {
            let mut z = Vec::new();
            let mut bins = 0;
            for (segment, per_knot) in segments.iter().zip(&samples) {
                if accept(segment) {
                    z.extend_from_slice(&per_knot[knot].0);
                    bins += per_knot[knot].1;
                }
            }
            (z, bins)
        };
        let mut beta2 = Vec::with_capacity(segments.len() * knots.len());
        for (segment, per_knot) in segments.iter().zip(&samples) {
            for (knot, (own, bins)) in per_knot.iter().enumerate() {
                let value = if *bins >= MIN_BINS {
                    variance(own)
                } else {
                    let same = |s: &Segment| ratio(s.cadence, segment.cadence) <= CADENCE_RATIO;
                    let (class, class_bins) = pooled(&same, knot);
                    if class_bins >= MIN_BINS {
                        variance(&class)
                    } else {
                        let (all, all_bins) = pooled(&|_| true, knot);
                        if all_bins >= MIN_BINS {
                            variance(&all)
                        } else {
                            1.0
                        }
                    }
                };
                beta2.push(value);
            }
        }
        Self { knots, beta2 }
    }

    /// Independent noise (beta = 1) for hand-built test samples.
    #[cfg(test)]
    pub fn white(segments: usize) -> Self {
        Self {
            knots: vec![1.0],
            beta2: vec![1.0; segments],
        }
    }

    /// Lower knot and the linear weight of the next knot for a window duration.
    pub fn interpolation(&self, duration: f64) -> (usize, f64) {
        let last = self.knots.len() - 1;
        if last == 0 || duration <= self.knots[0] {
            return (0, 0.0);
        }
        if duration >= self.knots[last] {
            return (last, 0.0);
        }
        let position =
            (duration / self.knots[0]).ln() / (self.knots[last] / self.knots[0]).ln() * last as f64;
        let lower = (position.floor() as usize).min(last - 1);
        (lower, position - lower as f64)
    }

    pub fn beta(&self, segment: usize, duration: f64) -> f64 {
        let (lower, fraction) = self.interpolation(duration);
        let row = &self.beta2[segment * self.knots.len()..(segment + 1) * self.knots.len()];
        let upper = row.get(lower + 1).copied().unwrap_or(row[lower]);
        ((1.0 - fraction) * row[lower] + fraction * upper).sqrt()
    }

    /// Per-point `w * beta^2`, knot-major: `[knot * points + point]`.
    pub fn point_weights(&self, weights: &[f64], segment_ids: &[usize]) -> Vec<f64> {
        let k = self.knots.len();
        (0..k)
            .flat_map(|knot| {
                weights
                    .iter()
                    .zip(segment_ids)
                    .map(move |(w, &s)| w * self.beta2[s * k + knot])
            })
            .collect()
    }
}

fn knots(min_duration: f64) -> Vec<f64> {
    if min_duration >= MAX_TRANSIT_DURATION {
        return vec![min_duration];
    }
    let step = (MAX_TRANSIT_DURATION / min_duration).ln() / (KNOTS - 1) as f64;
    (0..KNOTS)
        .map(|i| min_duration * (step * i as f64).exp())
        .collect()
}

/// beta^2 from normalized bin deviations: robust variance, clamped to [1, 100].
fn variance(z: &[f64]) -> f64 {
    robust_sigma(z).powi(2).clamp(1.0, MAX_BETA2)
}

/// Deviations of weighted bin means from the segment median, each divided by
/// its white-noise standard error, for `OFFSETS` shifted bin grids.
fn deviations(
    times: &[f64],
    values: &[f64],
    weights: &[f64],
    segment: &Segment,
    width: f64,
) -> (Vec<f64>, usize) {
    let range = segment.start..segment.end;
    // Zero-weight samples (masked transits) carry no information: they are
    // left out of the center, the bin means, and the bin point counts.
    let center = median(
        range
            .clone()
            .filter(|&i| weights[i] > 0.0)
            .map(|i| values[i]),
    );
    let origin = times[segment.start];
    let mut z = Vec::new();
    let mut independent = 0;
    for offset in 0..OFFSETS {
        let shift = offset as f64 / OFFSETS as f64;
        let mut current = None;
        let (mut sum, mut weight, mut count) = (0.0, 0.0, 0);
        let flush = |sum: f64, weight: f64, count: usize, z: &mut Vec<f64>| {
            if count >= MIN_BIN_POINTS && weight > 0.0 {
                z.push((sum / weight - center) * weight.sqrt());
                true
            } else {
                false
            }
        };
        for i in range.clone() {
            let bin = ((times[i] - origin) / width + shift).floor() as i64;
            if current != Some(bin) {
                if flush(sum, weight, count, &mut z) && offset == 0 {
                    independent += 1;
                }
                current = Some(bin);
                (sum, weight, count) = (0.0, 0.0, 0);
            }
            if weights[i] > 0.0 {
                sum += weights[i] * values[i];
                weight += weights[i];
                count += 1;
            }
        }
        if flush(sum, weight, count, &mut z) && offset == 0 {
            independent += 1;
        }
    }
    (z, independent)
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Noise(u64);

    impl Noise {
        fn gaussian(&mut self) -> f64 {
            let mut uniform = || {
                self.0 = self
                    .0
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                (((self.0 >> 11) as f64) + 0.5) / ((1_u64 << 53) as f64)
            };
            let (u, v) = (uniform(), uniform());
            (-2.0 * u.ln()).sqrt() * (std::f64::consts::TAU * v).cos()
        }
    }

    /// White noise plus an Ornstein-Uhlenbeck process with a 0.1-day timescale.
    fn series(count: usize, cadence: f64, white: f64, red: f64, seed: u64) -> Vec<f64> {
        let mut noise = Noise(seed);
        let decay = (-cadence / 0.1_f64).exp();
        let mut state = 0.0;
        (0..count)
            .map(|_| {
                state = state * decay + red * (1.0 - decay * decay).sqrt() * noise.gaussian();
                white * noise.gaussian() + state
            })
            .collect()
    }

    fn segment(start: usize, end: usize, cadence: f64) -> Segment {
        Segment {
            start,
            end,
            cadence,
        }
    }

    #[test]
    fn white_noise_has_beta_near_one_and_red_noise_grows_with_duration() {
        let cadence = 600.0 / 86400.0;
        let n = 4000;
        let times: Vec<f64> = (0..2 * n).map(|i| i as f64 * cadence).collect();
        let mut values = series(n, cadence, 400.0, 0.0, 1);
        values.extend(series(n, cadence, 400.0, 150.0, 2));
        let weights = vec![1.0 / 400.0_f64.powi(2); 2 * n];
        let segments = [segment(0, n, cadence), segment(n, 2 * n, cadence)];
        let red = RedNoise::estimate(&times, &values, &weights, &segments, 2.0 * cadence);
        assert_eq!(red.knots.len(), KNOTS);
        for knot in red.knots.iter().copied() {
            assert!(red.beta(0, knot) < 1.15, "white beta {}", red.beta(0, knot));
        }
        let short = red.beta(1, red.knots[0]);
        let long = red.beta(1, 2.0 / 24.0);
        assert!(short < 1.2, "{short}");
        // Expected near sqrt(1 + (0.75 * 150 / (400 / sqrt(12)))^2) ~ 1.35.
        assert!((1.2..1.6).contains(&long), "{long}");
    }

    #[test]
    fn interpolation_is_linear_in_log_duration_and_clamped() {
        let red = RedNoise {
            knots: knots(0.01),
            beta2: vec![1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 4.0],
        };
        assert_eq!(red.interpolation(0.001), (0, 0.0));
        assert_eq!(red.interpolation(1.0), (KNOTS - 1, 0.0));
        let between = (red.knots[6] * red.knots[7]).sqrt();
        let (lower, fraction) = red.interpolation(between);
        assert_eq!(lower, 6);
        assert!((fraction - 0.5).abs() < 1e-9);
        assert!((red.beta(0, between) - 2.5_f64.sqrt()).abs() < 1e-9);
        let single = RedNoise {
            knots: knots(0.5),
            beta2: vec![2.0],
        };
        assert_eq!(single.knots.len(), 1);
        assert!((single.beta(0, 0.1) - 2.0_f64.sqrt()).abs() < 1e-12);
    }

    #[test]
    fn short_segment_borrows_same_cadence_estimate() {
        let cadence = 600.0 / 86400.0;
        let mut values = series(4000, cadence, 400.0, 150.0, 3);
        values.extend(series(30, cadence, 400.0, 0.0, 4));
        let times: Vec<f64> = (0..4030)
            .map(|i| i as f64 * cadence + if i >= 4000 { 5.0 } else { 0.0 })
            .collect();
        let weights = vec![1.0 / 400.0_f64.powi(2); 4030];
        let segments = [segment(0, 4000, cadence), segment(4000, 4030, cadence)];
        let red = RedNoise::estimate(&times, &values, &weights, &segments, 2.0 * cadence);
        let last = KNOTS - 1;
        // 30 samples hold fewer than ten 0.35-day bins, so the class pool is used.
        let pooled = red.beta2[KNOTS + last];
        assert!((pooled / red.beta2[last] - 1.0).abs() < 0.2, "{pooled}");
        assert!(pooled > 1.0);
    }

    #[test]
    fn zero_weight_samples_are_equivalent_to_removed_samples() {
        let cadence = 600.0 / 86400.0;
        let n = 4000;
        let times: Vec<f64> = (0..n).map(|i| i as f64 * cadence).collect();
        let mut values = series(n, cadence, 400.0, 150.0, 5);
        // A deep masked "transit" stretch that must not influence beta.
        let masked: Vec<bool> = times.iter().map(|t| (10.0..10.3).contains(t)).collect();
        for (v, &m) in values.iter_mut().zip(&masked) {
            if m {
                *v -= 5000.0;
            }
        }
        let weights: Vec<f64> = masked
            .iter()
            .map(|&m| if m { 0.0 } else { 1.0 / 400.0_f64.powi(2) })
            .collect();
        let zeroed = RedNoise::estimate(
            &times,
            &values,
            &weights,
            &[segment(0, n, cadence)],
            2.0 * cadence,
        );
        let keep: Vec<usize> = (0..n).filter(|&i| !masked[i]).collect();
        let pick = |v: &[f64]| keep.iter().map(|&i| v[i]).collect::<Vec<f64>>();
        let removed = RedNoise::estimate(
            &pick(&times),
            &pick(&values),
            &pick(&weights),
            &[segment(0, keep.len(), cadence)],
            2.0 * cadence,
        );
        assert_eq!(zeroed.beta2, removed.beta2);
    }

    #[test]
    fn point_weights_scale_each_segment_by_its_beta_squared() {
        let red = RedNoise {
            knots: vec![0.1, 0.2],
            beta2: vec![1.0, 2.0, 3.0, 4.0],
        };
        let w = red.point_weights(&[1.0, 10.0], &[0, 1]);
        assert_eq!(w, vec![1.0, 30.0, 2.0, 40.0]);
    }
}

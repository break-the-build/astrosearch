use crate::{
    math::{median, robust_sigma},
    rednoise::RedNoise,
    segments::{self, Segment, CADENCE_RATIO},
    Point, TrendMethod,
};

/// Fewer neighbor differences than this borrow a pooled same-cadence estimate.
const MIN_NOISE_DIFFERENCES: usize = 20;

#[derive(Clone)]
pub(crate) struct Samples {
    pub origin: f64,
    pub times: Vec<f64>,
    pub transit: Vec<f64>,
    pub variability: Vec<f64>,
    pub transit_weights: Vec<f64>,
    pub variability_weights: Vec<f64>,
    /// Global median spacing; sets the shared search-grid resolution.
    pub cadence: f64,
    pub span: f64,
    /// Independent cadence/gap segments in time order.
    pub segments: Vec<Segment>,
    /// Segment index of each observation.
    pub segment_ids: Vec<usize>,
    /// Robust local scatter of each segment's transit series, in ppm.
    pub transit_noise: Vec<f64>,
    /// Duration-dependent correlated-noise factors of the transit series.
    pub red_noise: RedNoise,
}

#[cfg(test)]
pub(crate) fn prepare(
    points: &[Point],
    width: f64,
    method: TrendMethod,
) -> Result<Samples, String> {
    prepare_with_gap(points, width, width, method)
}

pub(crate) fn prepare_with_gap(
    points: &[Point],
    width: f64,
    gap_days: f64,
    method: TrendMethod,
) -> Result<Samples, String> {
    let cadence = median(points.windows(2).map(|p| p[1].time - p[0].time));
    let origin = points[0].time;
    let span = points.last().unwrap().time - origin;
    if !cadence.is_finite() || cadence <= 0.0 {
        return Err("observations must have a positive, finite cadence".into());
    }
    let absolute: Vec<f64> = points.iter().map(|p| p.time).collect();
    let segments = segments::split(&absolute, gap_days, cadence);
    let times: Vec<f64> = points.iter().map(|p| p.time - origin).collect();
    let long_width = (width.min(span) * 8.0).max(2.0).min(span);
    let short_trend = trend(points, &segments, width, method, None);
    let long_trend = trend(points, &segments, long_width, method, None);
    let mut transit = Vec::with_capacity(points.len());
    let mut variability = Vec::with_capacity(points.len());
    let mut transit_errors = Vec::with_capacity(points.len());
    let mut variability_errors = Vec::with_capacity(points.len());
    for (i, p) in points.iter().enumerate() {
        // The flux and its uncertainty use the same fitted baseline divisor.
        // Dividing errors by the instantaneous flux biases deep events' weights.
        transit.push((p.flux / short_trend[i] - 1.0) * 1e6);
        variability.push(p.flux / long_trend[i] - 1.0);
        transit_errors.push(p.error / short_trend[i] * 1e6);
        variability_errors.push(p.error / long_trend[i]);
    }
    if transit
        .iter()
        .chain(&variability)
        .chain(&transit_errors)
        .chain(&variability_errors)
        .any(|v| !v.is_finite())
    {
        return Err("flux or uncertainty dynamic range is too large to normalize".into());
    }
    let transit_noise = segment_noise(&transit, &times, &segments);
    let variability_noise = segment_noise(&variability, &times, &segments);
    let transit_weights = weights(&transit, &transit_errors, &segments, &transit_noise, 1e-6);
    let variability_weights = weights(
        &variability,
        &variability_errors,
        &segments,
        &variability_noise,
        1e-12,
    );
    let mut segment_ids = Vec::with_capacity(points.len());
    for (index, segment) in segments.iter().enumerate() {
        segment_ids.resize(segment.end, index);
    }
    // Transit windows start at two global cadences, the shortest searched duration.
    let red_noise =
        RedNoise::estimate(&times, &transit, &transit_weights, &segments, 2.0 * cadence);
    Ok(Samples {
        origin,
        times,
        transit,
        variability,
        transit_weights,
        variability_weights,
        cadence,
        span,
        segments,
        segment_ids,
        transit_noise,
        red_noise,
    })
}

/// Neighbor differences inside one segment, excluding gaps over three of its cadences.
fn differences(values: &[f64], times: &[f64], segment: &Segment) -> Vec<f64> {
    let range = segment.start..segment.end;
    values[range.clone()]
        .windows(2)
        .zip(times[range].windows(2))
        .filter(|(_, t)| t[1] - t[0] <= 3.0 * segment.cadence)
        .map(|(v, _)| v[1] - v[0])
        .collect()
}

/// Neighbor differences estimate local scatter without treating slow coherent
/// variability as independent noise. Each segment gets its own estimate, since
/// TESS sectors differ in cadence, camera, and scattered light. A segment with
/// too few differences borrows its cadence class, then the whole light curve.
/// Zero means no usable estimate.
fn segment_noise(values: &[f64], times: &[f64], segments: &[Segment]) -> Vec<f64> {
    let local: Vec<Vec<f64>> = segments
        .iter()
        .map(|s| differences(values, times, s))
        .collect();
    let pooled = |accept: &dyn Fn(&Segment) -> bool| -> Vec<f64> {
        segments
            .iter()
            .zip(&local)
            .filter(|(s, _)| accept(s))
            .flat_map(|(_, d)| d.iter().copied())
            .collect()
    };
    let scatter = |d: &[f64]| robust_sigma(d) / std::f64::consts::SQRT_2;
    let everywhere = scatter(&pooled(&|_| true));
    segments
        .iter()
        .zip(&local)
        .map(|(segment, own)| {
            if own.len() >= MIN_NOISE_DIFFERENCES {
                return scatter(own);
            }
            let class = pooled(&|s| segments::ratio(s.cadence, segment.cadence) <= CADENCE_RATIO);
            if class.len() >= MIN_NOISE_DIFFERENCES {
                scatter(&class)
            } else {
                everywhere
            }
        })
        .collect()
}

fn weights(
    values: &[f64],
    errors: &[f64],
    segments: &[Segment],
    noise: &[f64],
    minimum_sigma: f64,
) -> Vec<f64> {
    let global_fallback = robust_sigma(values);
    let mut output = Vec::with_capacity(values.len());
    for (segment, &local_noise) in segments.iter().zip(noise) {
        let fallback = if local_noise > 0.0 {
            local_noise
        } else {
            global_fallback
        };
        output.extend(errors[segment.start..segment.end].iter().map(|e| {
            let sigma = if *e > 0.0 {
                e.max(local_noise)
            } else {
                fallback
            };
            1.0 / sigma.max(minimum_sigma).powi(2)
        }));
    }
    output
}

/// Only occupied bins are stored, and independent segments never interpolate
/// across long gaps or cadence changes. Runtime and storage scale with
/// observations, not elapsed days.
fn trend(
    points: &[Point],
    segments: &[Segment],
    width: f64,
    method: TrendMethod,
    excluded: Option<&[bool]>,
) -> Vec<f64> {
    let mut output = Vec::with_capacity(points.len());
    for range in segments {
        let segment = &points[range.start..range.end];
        if method == TrendMethod::None {
            // Pre-flattened input: only remove each segment's median level.
            let level = median(
                segment
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| excluded.is_none_or(|mask| !mask[range.start + i]))
                    .map(|(_, p)| p.flux),
            );
            output.extend(std::iter::repeat_n(level, segment.len()));
            continue;
        }
        if method == TrendMethod::Biweight {
            let stride = (width / 4.0).max(range.cadence);
            let mut anchors = Vec::new();
            let mut center_time = segment[0].time;
            let last_time = segment.last().unwrap().time;
            while center_time < last_time {
                anchors.push(window_anchor(points, range, center_time, width, excluded));
                center_time += stride;
            }
            anchors.push(window_anchor(points, range, last_time, width, excluded));
            let mut anchor = 0;
            for point in segment {
                while anchor + 1 < anchors.len() && anchors[anchor + 1].0 <= point.time {
                    anchor += 1;
                }
                let value = if point.time <= anchors[0].0 || anchor + 1 == anchors.len() {
                    anchors[anchor].1
                } else {
                    let (ta, ya) = anchors[anchor];
                    let (tb, yb) = anchors[anchor + 1];
                    let fraction = (point.time - ta) / (tb - ta);
                    ya * (1.0 - fraction) + yb * fraction
                };
                output.push(value);
            }
            continue;
        }
        let mut anchors = Vec::new();
        let mut a = 0;
        while a < segment.len() {
            let mut b = a + 1;
            while b < segment.len() && segment[b].time - segment[a].time < width {
                b += 1;
            }
            let center = a + (b - a) / 2;
            let time = if (b - a).is_multiple_of(2) {
                0.5 * segment[center - 1].time + 0.5 * segment[center].time
            } else {
                segment[center].time
            };
            let values = segment[a..b]
                .iter()
                .enumerate()
                .filter(|(offset, _)| excluded.is_none_or(|mask| !mask[range.start + a + *offset]))
                .map(|(_, p)| p.flux);
            let baseline = match method {
                TrendMethod::Median | TrendMethod::None => median(values),
                TrendMethod::Biweight => tukey_biweight(values.collect()),
            };
            anchors.push((time, baseline));
            a = b;
        }
        let mut anchor = 0;
        for point in segment {
            while anchor + 1 < anchors.len() && anchors[anchor + 1].0 <= point.time {
                anchor += 1;
            }
            let value = if point.time <= anchors[0].0 || anchor + 1 == anchors.len() {
                anchors[anchor].1
            } else {
                let (ta, ya) = anchors[anchor];
                let (tb, yb) = anchors[anchor + 1];
                let fraction = (point.time - ta) / (tb - ta);
                // Convex form avoids overflow in subtraction for large fluxes.
                ya * (1.0 - fraction) + yb * fraction
            };
            output.push(value);
        }
    }
    output
}

pub(crate) fn detrend_masked(
    points: &[Point],
    samples: &Samples,
    width: f64,
    method: TrendMethod,
    excluded: &[bool],
) -> Vec<f64> {
    let baseline = trend(points, &samples.segments, width, method, Some(excluded));
    points
        .iter()
        .zip(baseline)
        .map(|(point, trend)| (point.flux / trend - 1.0) * 1e6)
        .collect()
}

fn window_anchor(
    points: &[Point],
    range: &Segment,
    center: f64,
    width: f64,
    excluded: Option<&[bool]>,
) -> (f64, f64) {
    let times =
        points[range.start..range.end].partition_point(|point| point.time < center - width / 2.0);
    let end =
        points[range.start..range.end].partition_point(|point| point.time <= center + width / 2.0);
    let left = range.start + times;
    let right = range.start + end;
    let values: Vec<f64> = (left..right)
        .filter(|&i| excluded.is_none_or(|mask| !mask[i]))
        .map(|i| points[i].flux)
        .collect();
    let value = if values.len() >= 5 {
        tukey_biweight(values)
    } else {
        let fallback: Vec<f64> = (left..right).map(|i| points[i].flux).collect();
        median(fallback)
    };
    (center, value)
}

/// Iteratively reweighted Tukey biweight location for a trend bin.
/// The redescending weights limit the effect of transits and isolated outliers.
fn tukey_biweight(values: Vec<f64>) -> f64 {
    if values.len() < 5 {
        return median(values);
    }
    let mut location = median(values.iter().copied());
    let scale = 1.4826 * median(values.iter().map(|value| (value - location).abs()));
    if !scale.is_finite() || scale <= f64::EPSILON * location.abs().max(1.0) {
        return location;
    }
    for _ in 0..5 {
        let mut weighted = 0.0;
        let mut weight_sum = 0.0;
        for value in &values {
            let u = (value - location) / (6.0 * scale);
            if u.abs() < 1.0 {
                let weight = (1.0 - u * u).powi(2);
                weighted += value * weight;
                weight_sum += weight;
            }
        }
        if weight_sum <= 0.0 {
            break;
        }
        let next = weighted / weight_sum;
        if (next - location).abs() <= 1e-10 * scale {
            location = next;
            break;
        }
        location = next;
    }
    location
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

    #[test]
    fn huge_gaps_do_not_allocate_empty_bins_or_mix_flux_offsets() {
        let points: Vec<Point> = (0..60)
            .map(|i| Point {
                time: if i < 30 {
                    i as f64 * 0.02
                } else {
                    1e9 + (i - 30) as f64 * 0.02
                },
                flux: if i < 30 { 1.0 } else { 100.0 },
                error: 0.0,
            })
            .collect();
        let prepared = prepare(&points, 0.5, TrendMethod::Median).unwrap();
        assert!(prepared.transit.iter().all(|v| v.abs() < 1e-6));
    }

    #[test]
    fn errors_are_normalized_by_trend_and_scatter_prevents_tiny_error_weights() {
        let mut points: Vec<Point> = (0..50)
            .map(|i| Point {
                time: i as f64 * 0.02,
                flux: 1.0,
                error: 0.001,
            })
            .collect();
        points[12].flux = 0.5;
        let prepared = prepare(&points, 0.5, TrendMethod::Median).unwrap();
        assert!((prepared.transit_weights[12] / prepared.transit_weights[11] - 1.0).abs() < 1e-9);
        let times: Vec<f64> = (0..50).map(|i| i as f64).collect();
        let values: Vec<f64> = (0..50).map(|i| (i as f64 * 2.17).sin()).collect();
        let segments = segments::split(&times, 2.0, 1.0);
        let noise = segment_noise(&values, &times, &segments);
        let w = weights(&values, &[1e-20; 50], &segments, &noise, 1e-6);
        assert!(w.iter().all(|v| *v < 100.0));
    }

    /// Equal noise per unit time: a 30-minute sample averages 15 two-minute ones.
    fn mixed_sectors(errors: bool) -> Vec<Point> {
        let mut noise = Noise(7);
        let mut points = Vec::new();
        for (start, seconds, count, sigma) in
            [(0.0, 1800.0, 600, 250e-6), (40.0, 120.0, 9000, 968e-6)]
        {
            for i in 0..count {
                points.push(Point {
                    time: start + i as f64 * seconds / 86400.0,
                    flux: 1.0 + sigma * noise.gaussian(),
                    error: if errors { sigma } else { 0.0 },
                });
            }
        }
        points
    }

    #[test]
    fn each_cadence_segment_gets_its_own_noise_and_weight() {
        for errors in [true, false] {
            let prepared = prepare(&mixed_sectors(errors), 0.5, TrendMethod::Median).unwrap();
            assert_eq!(prepared.segments.len(), 2);
            let [slow, fast] = [prepared.transit_noise[0], prepared.transit_noise[1]];
            assert!((slow / 250.0 - 1.0).abs() < 0.15, "slow noise {slow}");
            assert!((fast / 968.0 - 1.0).abs() < 0.15, "fast noise {fast}");
            // Sigma ratio ~3.9 means a variance-weight ratio near 15, not 1.
            let weight_ratio = prepared.transit_weights[0] / prepared.transit_weights[700];
            assert!(
                (11.0..20.0).contains(&weight_ratio),
                "{errors}: {weight_ratio}"
            );
        }
    }

    #[test]
    fn tiny_segment_borrows_same_cadence_noise_before_global_noise() {
        let mut points = mixed_sectors(false);
        // Five 30-minute samples, isolated by gaps, cannot estimate their own scatter.
        for i in 0..5 {
            points.push(Point {
                time: 80.0 + i as f64 * 1800.0 / 86400.0,
                flux: 1.0,
                error: 0.0,
            });
        }
        let prepared = prepare(&points, 0.5, TrendMethod::Median).unwrap();
        assert_eq!(prepared.segments.len(), 3);
        // The 30-minute pool (599 + 4 differences) matches sector 0, not the 2-minute sector.
        let [slow, fast, tiny] = [0, 1, 2].map(|i| prepared.transit_noise[i]);
        assert!((tiny / slow - 1.0).abs() < 0.02, "{tiny} vs {slow}");
        assert!(fast / tiny > 3.0, "{fast} vs {tiny}");
        assert!(prepared.segment_ids[..600].iter().all(|&s| s == 0));
        assert!(prepared.segment_ids[600..9600].iter().all(|&s| s == 1));
    }
}

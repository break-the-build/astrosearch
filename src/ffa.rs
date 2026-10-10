//! Fast folding algorithm (Staelin 1969) for the coarse stage of the
//! two-stage transit search.
//!
//! The series is placed on uniform slots of width `delta` (gaps are empty
//! slots). For a base period of `p` slots, the `m` rows of `p` slots are
//! combined by a radix-2 transform that yields folded profiles for all `m`
//! periods `p + s / (m - 1)` (s = 0..m) with `m p log2 m` additions instead of
//! `m * N`: neighbouring trial periods share partial sums. The slot width
//! grows with the shortest searched duration (half of it, in power-of-two
//! bands), following Ofir's (2014) duty-cycle-dependent binning, so long
//! periods are folded at coarse resolution. Profiles are scored with the
//! same white box statistic as `PhaseWorkspace::screen`; the best score of
//! each fine-grid cell ranks cells for exact refinement.

use crate::{
    grid::{maximum_duration, minimum_duration, DurationModel},
    rednoise::{RedNoise, MAX_TRANSIT_DURATION},
};
use std::{
    sync::atomic::{AtomicUsize, Ordering},
    thread,
};

/// Per-sample inputs of the box statistic; `times` start at zero.
pub(crate) struct FoldInputs<'a> {
    pub(crate) times: &'a [f64],
    pub(crate) sums: &'a [f64],
    pub(crate) background: &'a [f64],
    pub(crate) weights: &'a [f64],
}

/// What the coarse stage needs from the search configuration.
pub(crate) struct CoarsePlan<'a> {
    /// Full-resolution period grid; results are reported per cell of it.
    pub(crate) periods: &'a [f64],
    pub(crate) factor: usize,
    pub(crate) cadence: f64,
    pub(crate) model: DurationModel,
    pub(crate) durations: usize,
    pub(crate) threads: usize,
    /// One "segment" of transit-weight-averaged beta^2 per duration knot
    /// (all ones for white scoring), so windows are ranked with the same
    /// correlated statistic the exact search uses.
    pub(crate) red_noise: &'a RedNoise,
}

const QUANTITIES: usize = 4;
/// Per-slot clip of the mean residual, in units of that slot's noise.
const SLOT_CLIP_SIGMA: f64 = 5.0;

/// Uniform slots for one band: interleaved [sum, background, weight, count].
fn slot_series(inputs: &FoldInputs, delta: f64) -> Vec<f64> {
    let span = inputs.times.last().copied().unwrap_or(0.0);
    let slots = (span / delta).floor() as usize + 1;
    let mut out = vec![0.0; slots * QUANTITIES];
    for i in 0..inputs.times.len() {
        let w = inputs.weights[i];
        if w <= 0.0 {
            continue;
        }
        let s = ((inputs.times[i] / delta) as usize).min(slots - 1) * QUANTITIES;
        out[s] += inputs.sums[i];
        out[s + 1] += inputs.background[i];
        out[s + 2] += w;
        out[s + 3] += 1.0;
    }
    // A single deep artifact lands in some phase bin at every trial period;
    // the exact scan's repeated-event vetting rejects such windows, so cap
    // each slot's mean residual at SLOT_CLIP_SIGMA of its own noise. A real
    // transit contributes many moderate slots and is barely affected.
    for slot in out.chunks_exact_mut(QUANTITIES) {
        let w = slot[2];
        if w > 0.0 {
            let limit = SLOT_CLIP_SIGMA / w.sqrt();
            slot[0] = (slot[0] / w).clamp(-limit, limit) * w;
            slot[1] = (slot[1] / w).clamp(-limit, limit) * w;
        }
    }
    out
}

/// In-place radix-2 FFA over `m` rows (a power of two) of `p` interleaved
/// slots. Row `s` of the result is the fold at a drift of `s` slots across
/// all rows, i.e. period `p + s / (m - 1)` slots.
pub(crate) fn transform(data: &mut Vec<f64>, scratch: &mut Vec<f64>, m: usize, p: usize) {
    let row = p * QUANTITIES;
    scratch.resize(data.len(), 0.0);
    let mut n = 1;
    while n < m {
        for block in (0..m).step_by(2 * n) {
            for s in 0..2 * n {
                let (half, shift) = (s / 2, s.div_ceil(2));
                let a = (block + half) * row;
                let b = (block + n + half) * row;
                let out = (block + s) * row;
                for j in 0..p {
                    let source = ((j + shift) % p) * QUANTITIES;
                    let target = j * QUANTITIES;
                    for q in 0..QUANTITIES {
                        scratch[out + target + q] = data[a + target + q] + data[b + source + q];
                    }
                }
            }
        }
        std::mem::swap(data, scratch);
        n *= 2;
    }
}

/// Best box score of one folded profile over `widths` (slots, beta^2):
/// depth / sqrt(beta^2 / W_in + 1 / W_out), the exact search's correlated
/// statistic with one global beta^2 (beta^2 = 1 gives the white score).
fn best_score(
    profile: &[f64],
    p: usize,
    widths: &[(usize, f64)],
    prefix: &mut Vec<[f64; 4]>,
) -> f64 {
    prefix.clear();
    prefix.push([0.0; 4]);
    for i in 0..2 * p {
        let b = (i % p) * QUANTITIES;
        let last = prefix[i];
        prefix.push([
            last[0] + profile[b],
            last[1] + profile[b + 1],
            last[2] + profile[b + 2],
            last[3] + profile[b + 3],
        ]);
    }
    let total_background = prefix[p][1];
    let total_weight = prefix[p][2];
    let mut best = f64::NEG_INFINITY;
    for &(width, beta2) in widths {
        for start in 0..p {
            let (lo, hi) = (prefix[start], prefix[start + width]);
            let in_sum = hi[0] - lo[0];
            let in_weight = hi[2] - lo[2];
            let out_weight = total_weight - in_weight;
            if hi[3] - lo[3] < 4.0 || in_weight <= 0.0 || out_weight <= 0.0 {
                continue;
            }
            let out_background = total_background - (hi[1] - lo[1]);
            let numerator = out_background * in_weight - in_sum * out_weight;
            if numerator <= 0.0 {
                continue;
            }
            let depth = numerator / (in_weight * out_weight);
            let score = depth / (beta2.max(1.0) / in_weight + 1.0 / out_weight).sqrt();
            best = best.max(score);
        }
    }
    best
}

/// Window widths (in slots) and their beta^2 for the durations searched at `period`.
fn widths(period: f64, delta: f64, plan: &CoarsePlan, p: usize) -> Vec<(usize, f64)> {
    let shortest = minimum_duration(period, plan.cadence, plan.model);
    let longest = maximum_duration(period, plan.cadence, plan.model, MAX_TRANSIT_DURATION);
    let mut out: Vec<(usize, f64)> = (0..plan.durations)
        .map(|i| {
            let f = i as f64 / plan.durations.saturating_sub(1).max(1) as f64;
            let duration = shortest * (longest / shortest).powf(f);
            let width = ((duration / delta).round() as usize).clamp(1, (p / 3).max(1));
            (width, plan.red_noise.beta(0, width as f64 * delta).powi(2))
        })
        .collect();
    out.dedup_by_key(|w| w.0);
    out
}

/// One FFA job: a base period of `p` slots in a band of width `delta`.
struct Job {
    delta: f64,
    p: usize,
    band: usize,
}

/// Best coarse score of every fine-grid cell (NEG_INFINITY if unscored) and
/// the number of folded trial periods.
pub(crate) fn coarse_cell_scores(inputs: &FoldInputs, plan: &CoarsePlan) -> (Vec<f64>, usize) {
    let periods = plan.periods;
    let cells = periods.len().div_ceil(plan.factor);
    let (low, high) = (periods[0], periods[periods.len() - 1]);
    let base = plan
        .cadence
        .max(minimum_duration(low, plan.cadence, plan.model) / 2.0);
    // Power-of-two bands: band j uses slot width base * 2^j for periods whose
    // half minimum duration falls in [base 2^j, base 2^(j+1)).
    let mut bands: Vec<Vec<f64>> = Vec::new();
    let mut jobs = Vec::new();
    let mut delta = base;
    let mut start = low;
    while start < high {
        let mut end = start;
        while end < high && minimum_duration(end, plan.cadence, plan.model) / 2.0 < 2.0 * delta {
            end = (end * 1.01).min(high);
            if end >= high {
                break;
            }
        }
        let series = slot_series(inputs, delta);
        let slots = series.len() / QUANTITIES;
        let first = (start / delta).floor().max(2.0) as usize;
        let last = (end / delta).ceil() as usize;
        for p in first..=last {
            if slots >= 2 * p {
                jobs.push(Job {
                    delta,
                    p,
                    band: bands.len(),
                });
            }
        }
        bands.push(series);
        start = end;
        delta *= 2.0;
    }
    let workers = if plan.threads == 0 {
        thread::available_parallelism().map_or(1, usize::from)
    } else {
        plan.threads
    }
    .min(jobs.len().max(1));
    let next = AtomicUsize::new(0);
    let results: Vec<(Vec<f64>, usize)> = thread::scope(|scope| {
        let handles: Vec<_> = (0..workers)
            .map(|_| {
                let (jobs, bands, next) = (&jobs, &bands, &next);
                scope.spawn(move || {
                    let mut best = vec![f64::NEG_INFINITY; cells];
                    let (mut data, mut scratch, mut prefix) = (Vec::new(), Vec::new(), Vec::new());
                    let mut folded = 0;
                    loop {
                        let index = next.fetch_add(1, Ordering::Relaxed);
                        let Some(job) = jobs.get(index) else { break };
                        let series = &bands[job.band];
                        let slots = series.len() / QUANTITIES;
                        let rows = slots.div_ceil(job.p);
                        let m = rows.next_power_of_two();
                        data.clear();
                        data.resize(m * job.p * QUANTITIES, 0.0);
                        let used = slots * QUANTITIES;
                        data[..used].copy_from_slice(&series[..used]);
                        transform(&mut data, &mut scratch, m, job.p);
                        let row = job.p * QUANTITIES;
                        for s in 0..m {
                            let period = (job.p as f64 + s as f64 / (m - 1) as f64) * job.delta;
                            if period < low || period > high {
                                continue;
                            }
                            let w = widths(period, job.delta, plan, job.p);
                            let score =
                                best_score(&data[s * row..(s + 1) * row], job.p, &w, &mut prefix);
                            folded += 1;
                            let k = periods
                                .partition_point(|&q| q < period)
                                .min(periods.len() - 1);
                            let cell = k / plan.factor;
                            best[cell] = best[cell].max(score);
                        }
                    }
                    (best, folded)
                })
            })
            .collect();
        handles
            .into_iter()
            .map(|h| h.join().expect("ffa worker panicked"))
            .collect()
    });
    let mut best = vec![f64::NEG_INFINITY; cells];
    let mut folded = 0;
    for (partial, count) in results {
        folded += count;
        for (b, v) in best.iter_mut().zip(partial) {
            *b = b.max(v);
        }
    }
    (best, folded)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rows_with(values: &[f64], m: usize, p: usize) -> Vec<f64> {
        // Put one value per slot in the "sum" quantity; others zero.
        let mut data = vec![0.0; m * p * QUANTITIES];
        for (i, v) in values.iter().enumerate() {
            data[i * QUANTITIES] = *v;
        }
        data
    }

    #[test]
    fn transform_conserves_totals_and_aligns_drifting_pulses() {
        let (m, p) = (8, 16);
        // A pulse whose position drifts by exactly 1 slot per row: period p + 1/(m-1)*(m-1)/... ;
        // with drift s = m - 1 = 7 across the rows the fold must stack it in one bin.
        let mut values = vec![0.0; m * p];
        for r in 0..m {
            let shift = (r * 7) as f64 / (m - 1) as f64;
            values[r * p + (3 + shift.round() as usize) % p] = 1.0;
        }
        let mut data = rows_with(&values, m, p);
        let mut scratch = Vec::new();
        transform(&mut data, &mut scratch, m, p);
        let row = p * QUANTITIES;
        for s in 0..m {
            let total: f64 = (0..p).map(|j| data[s * row + j * QUANTITIES]).sum();
            assert!((total - m as f64).abs() < 1e-12, "row {s} total {total}");
        }
        let peak = (0..p)
            .map(|j| data[7 * row + j * QUANTITIES])
            .fold(0.0, f64::max);
        assert_eq!(peak, m as f64, "drift-7 row stacks every pulse");
        let unshifted = (0..p).map(|j| data[j * QUANTITIES]).fold(0.0, f64::max);
        assert!(unshifted < m as f64);
    }

    #[test]
    fn best_score_matches_direct_box_statistic() {
        let p = 20;
        let mut profile = vec![0.0; p * QUANTITIES];
        for j in 0..p {
            let dip = if (5..8).contains(&j) { -0.01 } else { 0.0 };
            profile[j * QUANTITIES] = dip * 10.0; // sum = w * residual, w = 10
            profile[j * QUANTITIES + 1] = dip * 10.0;
            profile[j * QUANTITIES + 2] = 10.0;
            profile[j * QUANTITIES + 3] = 5.0;
        }
        let score = best_score(&profile, p, &[(3, 1.0)], &mut Vec::new());
        // depth 0.01 over in-weight 30 and out-weight 170.
        let expected = 0.01 / (1.0 / 30.0 + 1.0 / 170.0f64).sqrt();
        assert!((score - expected).abs() < 1e-9, "{score} vs {expected}");
        // Correlated noise (beta^2 = 4) inflates the in-transit variance only.
        let red = best_score(&profile, p, &[(3, 4.0)], &mut Vec::new());
        let expected_red = 0.01 / (4.0 / 30.0 + 1.0 / 170.0f64).sqrt();
        assert!((red - expected_red).abs() < 1e-9, "{red} vs {expected_red}");
    }
}

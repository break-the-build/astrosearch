//! Per-worker phase-folding buffers for one trial period of the box search.

use crate::epochs::CycleRange;

#[derive(Default)]
pub(crate) struct PhaseWorkspace {
    pub(crate) sums: Vec<f64>,
    pub(crate) background: Vec<f64>,
    pub(crate) weights: Vec<f64>,
    pub(crate) counts: Vec<usize>,
    /// Per duration knot, the phase histogram of `w * beta^2` (knot-major).
    pub(crate) red: Vec<f64>,
    pub(crate) sum_prefix: Vec<f64>,
    pub(crate) background_prefix: Vec<f64>,
    pub(crate) weight_prefix: Vec<f64>,
    pub(crate) count_prefix: Vec<usize>,
    pub(crate) red_prefix: Vec<f64>,
    /// Which knots' histograms are built for the current trial period.
    red_built: Vec<bool>,
    pub(crate) point_bins: Vec<usize>,
    pub(crate) point_cycles: Vec<i64>,
    pub(crate) point_sum_prefix: Vec<f64>,
    pub(crate) point_weight_prefix: Vec<f64>,
    /// Running count of positive-weight (unmasked) samples.
    pub(crate) point_count_prefix: Vec<usize>,
    pub(crate) cycles: Vec<CycleRange>,
    /// Correlated-noise histograms and per-event bookkeeping are needed only
    /// for windows whose white score clears the threshold, which most trial
    /// periods never have, so both are built on first use.
    red_ready: bool,
    epochs_ready: bool,
}

impl PhaseWorkspace {
    pub(crate) fn reset(&mut self, bins: usize, observations: usize) {
        for values in [&mut self.sums, &mut self.background, &mut self.weights] {
            values.resize(bins, 0.0);
            values.fill(0.0);
        }
        self.counts.resize(bins, 0);
        self.counts.fill(0);
        for values in [
            &mut self.sum_prefix,
            &mut self.background_prefix,
            &mut self.weight_prefix,
        ] {
            values.resize(2 * bins + 1, 0.0);
            values[0] = 0.0;
        }
        self.count_prefix.resize(2 * bins + 1, 0);
        self.count_prefix[0] = 0;
        self.point_bins.resize(observations, 0);
        self.point_cycles.resize(observations, 0);
        self.red_ready = false;
        self.epochs_ready = false;
    }

    /// Circular prefix sums over two laps of the phase bins, so any window
    /// `start..start + width` with `start < bins` is a single difference.
    pub(crate) fn prefix(&mut self, bins: usize) {
        for i in 0..2 * bins {
            // Two laps without a per-element modulo.
            let b = if i < bins { i } else { i - bins };
            self.sum_prefix[i + 1] = self.sum_prefix[i] + self.sums[b];
            self.background_prefix[i + 1] = self.background_prefix[i] + self.background[b];
            self.weight_prefix[i + 1] = self.weight_prefix[i] + self.weights[b];
            self.count_prefix[i + 1] = self.count_prefix[i] + self.counts[b];
        }
    }

    /// Window starts whose white score could reach `floor`, in phase order.
    /// A tight, division-free pass using `n = depth Win Wout`: white >= floor
    /// exactly when n > 0 and n^2 >= floor^2 Win Wout (Win + Wout). Windows are
    /// dropped only when they fail by a 1e-9 relative margin; survivors are
    /// re-scored exactly by the caller, so screening never changes a result.
    pub(crate) fn screen(&self, bins: usize, width: usize, floor: f64, survivors: &mut Vec<usize>) {
        survivors.clear();
        let total_weight = self.weight_prefix[bins];
        let total_background = self.background_prefix[bins];
        let limit = floor * floor * (1.0 - 1e-9);
        let count = &self.count_prefix;
        let weight = &self.weight_prefix;
        let sum = &self.sum_prefix;
        let background = &self.background_prefix;
        for start in 0..bins {
            let end = start + width;
            let in_weight = weight[end] - weight[start];
            let out_weight = total_weight - in_weight;
            let in_sum = sum[end] - sum[start];
            let out_background = total_background - (background[end] - background[start]);
            let numerator = out_background * in_weight - in_sum * out_weight;
            // Non-short-circuit `&`: every term is cheap and branch-free
            // evaluation lets the compiler vectorize the loop body.
            let possible = (count[end] - count[start] >= 4)
                & (in_weight > 0.0)
                & (out_weight > 0.0)
                & (numerator > 0.0)
                & (numerator * numerator
                    >= limit * in_weight * out_weight * (in_weight + out_weight));
            if possible {
                survivors.push(start);
            }
        }
    }

    /// Phase histogram of `w * beta^2` for duration knot `knot` and its
    /// circular prefix sums, built on first use within a trial period
    /// (`point_weights` is knot-major with `knots` rows). A window reads only
    /// its own knot and the next one, so most knots are never built.
    pub(crate) fn ensure_red_knot(
        &mut self,
        bins: usize,
        knots: usize,
        knot: usize,
        point_weights: &[f64],
    ) {
        if !self.red_ready {
            self.red.resize(knots * bins, 0.0);
            self.red_prefix.resize(knots * (2 * bins + 1), 0.0);
            self.red_built.clear();
            self.red_built.resize(knots, false);
            self.red_ready = true;
        }
        if self.red_built[knot] {
            return;
        }
        let observations = self.point_bins.len();
        let histogram = &mut self.red[knot * bins..(knot + 1) * bins];
        histogram.fill(0.0);
        // Each bin accumulates in sample order, as in a sample-major loop.
        let weights = &point_weights[knot * observations..(knot + 1) * observations];
        for (&b, &w) in self.point_bins.iter().zip(weights) {
            histogram[b] += w;
        }
        let stride = 2 * bins + 1;
        let prefix = &mut self.red_prefix[knot * stride..(knot + 1) * stride];
        prefix[0] = 0.0;
        for i in 0..2 * bins {
            prefix[i + 1] = prefix[i] + histogram[if i < bins { i } else { i - bins }];
        }
        self.red_built[knot] = true;
    }

    /// Time-ordered prefix sums and per-cycle index ranges for event statistics.
    pub(crate) fn ensure_epochs(&mut self, sums: &[f64], weights: &[f64]) {
        if self.epochs_ready {
            return;
        }
        let observations = self.point_bins.len();
        self.point_sum_prefix.resize(observations + 1, 0.0);
        self.point_weight_prefix.resize(observations + 1, 0.0);
        self.point_count_prefix.resize(observations + 1, 0);
        self.point_sum_prefix[0] = 0.0;
        self.point_weight_prefix[0] = 0.0;
        self.point_count_prefix[0] = 0;
        self.cycles.clear();
        let (mut point_sum, mut point_weight, mut point_count) = (0.0, 0.0, 0);
        for i in 0..observations {
            point_sum += sums[i];
            point_weight += weights[i];
            point_count += usize::from(weights[i] > 0.0);
            self.point_sum_prefix[i + 1] = point_sum;
            self.point_weight_prefix[i + 1] = point_weight;
            self.point_count_prefix[i + 1] = point_count;
            let cycle = self.point_cycles[i];
            if let Some(range) = self.cycles.last_mut().filter(|r| r.cycle == cycle) {
                range.end = i + 1;
            } else {
                self.cycles.push(CycleRange {
                    cycle,
                    start: i,
                    end: i + 1,
                });
            }
        }
        self.epochs_ready = true;
    }

    /// Sum of `w * beta^2` at one knot over circular bins `start..end`.
    pub(crate) fn red_sum(&self, bins: usize, knot: usize, start: usize, end: usize) -> f64 {
        let offset = knot * (2 * bins + 1);
        self.red_prefix[offset + end] - self.red_prefix[offset + start]
    }
}

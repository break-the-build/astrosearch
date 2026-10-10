//! Weighted, phase-folded box search for periodic dips.

use crate::{
    epochs::{repeated_evidence, Epoch, EpochIndex},
    grid::{
        adaptive_period_grid, maximum_duration, minimum_duration, minimum_orbital_period,
        DurationModel, TrialBudget,
    },
    grouping::{collapse_periodic, collapse_transit_events, event_support_overlap},
    math::{median, robust_sigma},
    parallel::collect_trials,
    phase::PhaseWorkspace,
    preprocess::Samples,
    rednoise::{RedNoise, MAX_TRANSIT_DURATION},
    Candidate, CandidateKind, Config, EventEvidence, Scoring,
};

pub(crate) struct TransitSearch {
    /// Transit candidates, strongest first, after the configured grouping rule.
    pub candidates: Vec<Candidate>,
    /// Best accepted window of every trial period, when requested.
    pub periodogram: Vec<Candidate>,
    pub trials: usize,
}

const MAX_CANDIDATES_PER_PERIOD: usize = 3;

/// Iteratively remove the strongest transit windows and rerun the BLS scan to
/// expose additional planets hidden below the first signal and its aliases.
///
/// Each pass contributes its strongest candidate. The remaining output slots
/// are filled from the last pass's ranked list, which was searched with every
/// earlier winner masked, so `transit_candidates` still bounds the output
/// rather than the number of passes.
pub(crate) fn iterative_periodic_dips(
    samples: &Samples,
    config: Config,
    warnings: &mut Vec<String>,
) -> TransitSearch {
    let mut masked = samples.clone();
    let mut candidates = Vec::new();
    let mut remaining = Vec::new();
    let mut periodogram = Vec::new();
    let mut trials = 0;
    for pass in 0..config.transit_search_passes.min(config.transit_candidates) {
        let search = periodic_dips(&masked, config, warnings);
        trials += search.trials;
        if pass == 0 {
            periodogram = search.periodogram;
        }
        let mut ranked = search.candidates.into_iter();
        let Some(candidate) = ranked.next() else {
            remaining.clear();
            break;
        };
        if candidate.event_evidence.is_empty() {
            // Nothing to mask: this pass's whole list is the final remainder.
            remaining = std::iter::once(candidate).chain(ranked).collect();
            break;
        }
        remaining = ranked.collect();
        if config.stop_passes_below_detection && candidate.score < config.detection_threshold {
            // Jackknife <= score, so nothing in later passes could reach the
            // detection tier; keep this pass's list and stop.
            candidates.push(candidate);
            break;
        }
        // Masked samples are excluded (weight zero) rather than replaced by a
        // baseline value, so they neither dilute later windows nor count as data.
        let mut newly_masked = 0usize;
        for (index, &time) in masked.times.iter().enumerate() {
            let absolute_time = samples.origin + time;
            let in_event = candidate.event_evidence.iter().any(|event| {
                (absolute_time - event.epoch).abs() <= candidate.duration / 2.0 + samples.cadence
            });
            if in_event && masked.transit_weights[index] > 0.0 {
                masked.transit_weights[index] = 0.0;
                newly_masked += 1;
            }
        }
        candidates.push(candidate);
        if newly_masked == 0 {
            break;
        }
        // The removed transit no longer inflates the correlated-noise estimate.
        masked.red_noise = RedNoise::estimate(
            &masked.times,
            &masked.transit,
            &masked.transit_weights,
            &masked.segments,
            2.0 * masked.cadence,
        );
    }
    candidates.extend(remaining);
    let mut candidates = match config.candidate_grouping {
        crate::CandidateGrouping::Ephemeris => {
            collapse_periodic(candidates, samples.span, config.transit_candidates)
        }
        crate::CandidateGrouping::Events => {
            collapse_transit_events(candidates, config.transit_candidates)
        }
    };
    // A periodic candidate needs three or more observed transits; with two,
    // the period is only known up to an integer fraction of their separation.
    for candidate in candidates.iter_mut().filter(|c| c.n_events < 3) {
        candidate.kind = CandidateKind::DuoTransit;
        candidate
            .note
            .push_str("; two observed transits: period ambiguous, enumerate aliases");
    }
    TransitSearch {
        candidates,
        periodogram,
        trials,
    }
}

/// A distinct phase-window width and its correlated-noise interpolation.
#[derive(Clone, Copy)]
struct Window {
    width: usize,
    knot: usize,
    fraction: f64,
}

/// Correlated scoring inflates the in-window variance from 1/W to
/// sum(w * beta^2)/W^2, using each sample's segment beta at the window's
/// duration. Because beta >= 1, it never exceeds the white score, so the white
/// score safely prunes windows before the correlated one decides selection.
pub(crate) fn periodic_dips(
    samples: &Samples,
    config: Config,
    warnings: &mut Vec<String>,
) -> TransitSearch {
    let min_duration = minimum_duration(
        config.min_period,
        samples.cadence,
        DurationModel::from_config(&config),
    );
    let model = DurationModel::from_config(&config);
    let mut min_period = config.min_period.max(3.0 * min_duration);
    if model.physical {
        min_period = min_period.max(minimum_orbital_period(model));
    }
    // Two events can be observed at periods longer than half the baseline.
    // Actual covered cycles are checked below, rather than vetoed by period.
    let max_period = config.max_period.min(samples.span);
    if max_period <= min_period {
        return TransitSearch {
            candidates: Vec::new(),
            periodogram: Vec::new(),
            trials: 0,
        };
    }
    let periods = adaptive_period_grid(
        min_period,
        max_period,
        samples.span,
        samples.cadence,
        DurationModel::from_config(&config),
        TrialBudget::from_config(&config),
        warnings,
    );
    let trials = periods.len();
    // Masked samples (weight zero) from earlier search passes are ignored.
    let unmasked: Vec<f64> = samples
        .transit
        .iter()
        .zip(&samples.transit_weights)
        .filter(|(_, &w)| w > 0.0)
        .map(|(&v, _)| v)
        .collect();
    let baseline = median(unmasked.iter().copied());
    let clip = (5.0 * robust_sigma(&unmasked)).max(1e-6);
    // Reuse weighted residuals and phase/prefix allocations across periods.
    let residuals: Vec<f64> = samples.transit.iter().map(|v| v - baseline).collect();
    let sums: Vec<f64> = residuals
        .iter()
        .zip(&samples.transit_weights)
        .map(|(y, w)| y * w)
        .collect();
    let background: Vec<f64> = residuals
        .iter()
        .zip(&samples.transit_weights)
        .map(|(y, w)| y.clamp(-clip, clip) * w)
        .collect();
    let correlated = config.scoring == Scoring::Correlated;
    let red_weights = if correlated {
        samples
            .red_noise
            .point_weights(&samples.transit_weights, &samples.segment_ids)
    } else {
        Vec::new()
    };
    let context = TrialContext {
        samples,
        config,
        periods: &periods,
        sums: &sums,
        background: &background,
        red_weights: &red_weights,
        correlated,
    };
    let (per_trial, evaluated) = if config.two_stage_factor >= 2 {
        two_stage_scan(&context, warnings)
    } else {
        let all = collect_trials(
            trials,
            config.threads,
            TrialScratch::default,
            |scratch, k| context.scan(k, scratch),
        );
        (all, trials)
    };
    let mut periodogram = Vec::new();
    if config.keep_periodogram {
        let mut last_period = None;
        for candidate in &per_trial {
            if last_period != Some(candidate.period.to_bits()) {
                periodogram.push(candidate.clone());
                last_period = Some(candidate.period.to_bits());
            }
        }
    }
    let candidates = match config.candidate_grouping {
        crate::CandidateGrouping::Ephemeris => {
            collapse_periodic(per_trial, samples.span, config.transit_candidates)
        }
        crate::CandidateGrouping::Events => {
            collapse_transit_events(per_trial, config.transit_candidates)
        }
    };
    TransitSearch {
        candidates,
        periodogram,
        trials: evaluated,
    }
}

/// Score floor of the coarse stage, as a fraction of `min_transit_score`:
/// a coarse period is up to ~factor/4 minimum durations off, which smears a
/// transit over more phase and lowers its score.
const COARSE_FLOOR_FRACTION: f64 = 0.6;

/// Transit-weight-averaged beta^2 per duration knot as a one-segment model
/// (white when correlated scoring is off).
fn global_red_noise(samples: &Samples, correlated: bool) -> RedNoise {
    let red = &samples.red_noise;
    let knots = red.knots.len();
    if !correlated || red.beta2.len() < knots {
        return RedNoise::white(1);
    }
    let segments = red.beta2.len() / knots;
    let mut weight = vec![0.0; segments];
    for (&w, &s) in samples.transit_weights.iter().zip(&samples.segment_ids) {
        weight[s.min(segments - 1)] += w;
    }
    let total: f64 = weight.iter().sum();
    if total <= 0.0 {
        return RedNoise::white(1);
    }
    let beta2 = (0..knots)
        .map(|k| {
            (0..segments)
                .map(|s| weight[s] * red.beta2[s * knots + k])
                .sum::<f64>()
                / total
        })
        .collect();
    RedNoise {
        knots: red.knots.clone(),
        beta2,
    }
}

/// Coarse-to-fine search over the full-resolution grid in `context.periods`.
///
/// Stage 1 evaluates the centre of every cell of `two_stage_factor` adjacent
/// fine periods with a lowered score floor. Stage 2 evaluates every fine
/// period of the `two_stage_keep` best-scoring cells and their neighbours with
/// the normal settings; only stage-2 results are returned. Selection is by
/// score with index tie-breaks, so results are deterministic and independent
/// of the thread count. Returns the candidates and the number of trial
/// periods evaluated.
fn two_stage_scan(context: &TrialContext, warnings: &mut Vec<String>) -> (Vec<Candidate>, usize) {
    let config = context.config;
    let fine = context.periods.len();
    let factor = config.two_stage_factor;
    let cells = fine.div_ceil(factor);
    let centre = |cell: usize| (cell * factor + factor / 2).min(fine - 1);
    let floor = config.min_transit_score * COARSE_FLOOR_FRACTION;
    let (best, coarse_trials) = if config.coarse_ffa {
        let inputs = crate::ffa::FoldInputs {
            times: &context.samples.times,
            sums: context.sums,
            background: context.background,
            weights: &context.samples.transit_weights,
        };
        let red_noise = global_red_noise(context.samples, context.correlated);
        let plan = crate::ffa::CoarsePlan {
            periods: context.periods,
            factor,
            cadence: context.samples.cadence,
            model: DurationModel::from_config(&config),
            durations: config.durations,
            threads: config.threads,
            red_noise: &red_noise,
        };
        crate::ffa::coarse_cell_scores(&inputs, &plan)
    } else {
        let coarse_context = TrialContext {
            config: Config {
                min_transit_score: floor,
                ..config
            },
            ..*context
        };
        let coarse = collect_trials(
            cells,
            config.threads,
            TrialScratch::default,
            |scratch, cell| coarse_context.scan(centre(cell), scratch),
        );
        let index_of: std::collections::HashMap<u64, usize> = (0..cells)
            .map(|cell| (context.periods[centre(cell)].to_bits(), cell))
            .collect();
        let mut best = vec![f64::NEG_INFINITY; cells];
        for candidate in &coarse {
            if let Some(&cell) = index_of.get(&candidate.period.to_bits()) {
                best[cell] = best[cell].max(candidate.score);
            }
        }
        (best, cells)
    };
    let mut ranked: Vec<usize> = (0..cells).filter(|&c| best[c] >= floor).collect();
    ranked.sort_by(|&a, &b| best[b].total_cmp(&best[a]).then(a.cmp(&b)));
    ranked.truncate(config.two_stage_keep);
    let mut refine = vec![false; fine];
    for &cell in &ranked {
        let start = cell.saturating_sub(1) * factor;
        let end = ((cell + 2) * factor).min(fine);
        refine[start..end].fill(true);
    }
    let indices: Vec<usize> = (0..fine).filter(|&k| refine[k]).collect();
    let refined = collect_trials(
        indices.len(),
        config.threads,
        TrialScratch::default,
        |scratch, i| context.scan(indices[i], scratch),
    );
    warnings.push(format!(
        "two-stage search ({}): {coarse_trials} coarse trials, {} cells kept, {} of {fine} fine trials refined",
        if config.coarse_ffa { "fast folding" } else { "grid" },
        ranked.len(),
        indices.len()
    ));
    (refined, coarse_trials + indices.len())
}

/// Read-only inputs shared by every trial period of one search pass.
#[derive(Clone, Copy)]
struct TrialContext<'a> {
    samples: &'a Samples,
    config: Config,
    periods: &'a [f64],
    /// Per-sample `w * residual` and `w * clipped residual`.
    sums: &'a [f64],
    background: &'a [f64],
    /// Per-sample `w * beta^2` at each duration knot (empty for white scoring).
    red_weights: &'a [f64],
    correlated: bool,
}

/// Per-worker buffers reused across trial periods.
#[derive(Default)]
struct TrialScratch {
    workspace: PhaseWorkspace,
    epochs: Vec<Epoch>,
    wrapped: Vec<Epoch>,
    merged: Vec<Epoch>,
    windows: Vec<Window>,
    survivors: Vec<usize>,
}

/// Phase-bin geometry of one trial period.
#[derive(Clone, Copy)]
struct Trial {
    period: f64,
    bins: usize,
    red_knots: usize,
}

impl TrialContext<'_> {
    /// Best accepted windows of trial period `k`: one, or a short list of
    /// distinct-event windows when grouping by events.
    fn scan(&self, k: usize, scratch: &mut TrialScratch) -> Vec<Candidate> {
        let config = self.config;
        let period = self.periods[k];
        let candidate_limit = match config.candidate_grouping {
            crate::CandidateGrouping::Ephemeris => 1,
            crate::CandidateGrouping::Events => MAX_CANDIDATES_PER_PERIOD,
        };
        let period_min_duration = minimum_duration(
            period,
            self.samples.cadence,
            DurationModel::from_config(&config),
        );
        let bins = ((2.0 * period / period_min_duration).ceil() as usize)
            .clamp(128, 16_384)
            .next_power_of_two();
        let TrialScratch {
            workspace,
            epochs,
            wrapped,
            merged,
            windows,
            survivors,
        } = scratch;
        self.plan_windows(period, bins, period_min_duration, windows);
        let red_knots = windows
            .iter()
            .filter(|_| self.correlated)
            .map(|w| w.knot + 1 + usize::from(w.fraction > 0.0))
            .max()
            .unwrap_or(0);
        let trial = Trial {
            period,
            bins,
            red_knots,
        };
        self.fold(trial, workspace);
        let mut best: Vec<Candidate> = Vec::with_capacity(MAX_CANDIDATES_PER_PERIOD);
        for window in windows.iter() {
            // The fixed threshold never exceeds the evolving floor below,
            // so the screen keeps every window the exact logic could accept.
            workspace.screen(bins, window.width, config.min_transit_score, survivors);
            for &start in survivors.iter() {
                let weakest = (best.len() == candidate_limit)
                    .then(|| best.last().map(|candidate| candidate.score).unwrap_or(0.0));
                let buffers = (&mut *epochs, &mut *wrapped, &mut *merged);
                if let Some(candidate) =
                    self.evaluate(trial, *window, start, weakest, workspace, buffers)
                {
                    keep_best(&mut best, candidate, candidate_limit);
                }
            }
        }
        // Keep only the strongest windows per trial period; event mode
        // retains a small shortlist so distinct observed events can survive.
        best
    }

    /// Distinct phase-window widths for the searched durations.
    fn plan_windows(
        &self,
        period: f64,
        bins: usize,
        period_min_duration: f64,
        windows: &mut Vec<Window>,
    ) {
        let config = self.config;
        let max_duration = maximum_duration(
            period,
            self.samples.cadence,
            DurationModel::from_config(&config),
            MAX_TRANSIT_DURATION,
        )
        .max(period_min_duration);
        windows.clear();
        for i in 0..config.durations {
            let fraction = i as f64 / config.durations.saturating_sub(1).max(1) as f64;
            let duration =
                period_min_duration * (max_duration / period_min_duration).powf(fraction);
            let width = ((duration / period * bins as f64).round() as usize).clamp(1, bins / 3);
            // Different durations often round to the same phase window.
            if windows.last().is_some_and(|w| w.width == width) {
                continue;
            }
            let (knot, fraction) = if self.correlated {
                self.samples
                    .red_noise
                    .interpolation(width as f64 / bins as f64 * period)
            } else {
                (0, 0.0)
            };
            windows.push(Window {
                width,
                knot,
                fraction,
            });
        }
    }

    /// Fold every sample into phase bins and build the circular prefix sums.
    fn fold(&self, trial: Trial, workspace: &mut PhaseWorkspace) {
        let (samples, period, bins) = (self.samples, trial.period, trial.bins);
        workspace.reset(bins, samples.times.len());
        for (i, time) in samples.times.iter().enumerate() {
            let position = time / period;
            let cycle = position.floor();
            // Times are measured from the first sample, so position >= 0,
            // where position - floor(position) is exact and bitwise equal
            // to rem_euclid(1.0) without its fmod call.
            let phase = if position >= 0.0 {
                position - cycle
            } else {
                position.rem_euclid(1.0)
            };
            let b = ((phase * bins as f64) as usize).min(bins - 1);
            workspace.point_bins[i] = b;
            workspace.point_cycles[i] = cycle as i64;
            workspace.sums[b] += self.sums[i];
            workspace.background[b] += self.background[i];
            workspace.weights[b] += samples.transit_weights[i];
            workspace.counts[b] += usize::from(samples.transit_weights[i] > 0.0);
        }
        workspace.prefix(bins);
    }

    /// Score one window and, if it beats `weakest` and passes event vetting,
    /// return it as a candidate.
    fn evaluate(
        &self,
        trial: Trial,
        window: Window,
        start: usize,
        weakest: Option<f64>,
        workspace: &mut PhaseWorkspace,
        (epochs, wrapped, merged): (&mut Vec<Epoch>, &mut Vec<Epoch>, &mut Vec<Epoch>),
    ) -> Option<Candidate> {
        let (config, samples) = (self.config, self.samples);
        let (period, bins, width) = (trial.period, trial.bins, window.width);
        let end = start + width;
        let count = workspace.count_prefix[end] - workspace.count_prefix[start];
        let in_weight = workspace.weight_prefix[end] - workspace.weight_prefix[start];
        let out_weight = workspace.weight_prefix[bins] - in_weight;
        if count < 4 || in_weight <= 0.0 || out_weight <= 0.0 {
            return None;
        }
        let in_sum = workspace.sum_prefix[end] - workspace.sum_prefix[start];
        let out_background = workspace.background_prefix[bins]
            - (workspace.background_prefix[end] - workspace.background_prefix[start]);
        // Division-free screen: white >= floor exactly when
        // n > 0 and n^2 >= floor^2 Win Wout (Win + Wout), with
        // n = depth Win Wout. Only windows failing by a 1e-9 relative
        // margin are skipped, so every window that could pass is
        // still scored by the exact expression below.
        let floor = weakest.map_or(config.min_transit_score, |w| {
            w.max(config.min_transit_score)
        });
        let numerator = out_background * in_weight - in_sum * out_weight;
        if numerator.is_nan()
            || numerator <= 0.0
            || numerator * numerator
                < floor * floor * in_weight * out_weight * (in_weight + out_weight) * (1.0 - 1e-9)
        {
            return None;
        }
        let in_mean = in_sum / in_weight;
        let out_mean = out_background / out_weight;
        let depth = out_mean - in_mean;
        // Robust two-level SNR proxy: clipping the background means it
        // is not the exact Gaussian likelihood-ratio statistic.
        let white = depth / (1.0 / in_weight + 1.0 / out_weight).sqrt();
        let beaten = |score: f64| {
            !score.is_finite()
                || score < config.min_transit_score
                || weakest.is_some_and(|minimum| score <= minimum)
        };
        if beaten(white) {
            return None;
        }
        let (score, noise_factor) = if self.correlated {
            workspace.ensure_red_knot(bins, trial.red_knots, window.knot, self.red_weights);
            let mut red_weight = workspace.red_sum(bins, window.knot, start, end);
            if window.fraction > 0.0 {
                let upper_knot = window.knot + 1;
                workspace.ensure_red_knot(bins, trial.red_knots, upper_knot, self.red_weights);
                let upper = workspace.red_sum(bins, upper_knot, start, end);
                red_weight = (1.0 - window.fraction) * red_weight + window.fraction * upper;
            }
            let in_variance = (red_weight / (in_weight * in_weight)).max(1.0 / in_weight);
            let score = depth / (in_variance + 1.0 / out_weight).sqrt();
            if beaten(score) {
                return None;
            }
            (score, Some(white / score))
        } else {
            (white, None)
        };
        workspace.ensure_epochs(self.sums, &samples.transit_weights);
        let epoch_index = EpochIndex {
            cycles: &workspace.cycles,
            point_bins: &workspace.point_bins,
            point_sum_prefix: &workspace.point_sum_prefix,
            point_weight_prefix: &workspace.point_weight_prefix,
            point_count_prefix: &workspace.point_count_prefix,
        };
        epoch_index.collect_window_epochs(start, width, bins, epochs, wrapped, merged);
        let (n_events, odd_even_sigma) = repeated_evidence(epochs, out_mean, out_weight)?;
        let phase_center = ((start as f64 + width as f64 / 2.0) / bins as f64).rem_euclid(1.0);
        let event_evidence = epochs
            .iter()
            .filter(|event| event.points >= 2 && event.weight > 0.0)
            .map(|event| {
                let depth = out_mean - event.sum / event.weight;
                let significance = depth / (1.0 / event.weight + 1.0 / out_weight).sqrt();
                EventEvidence {
                    epoch: samples.origin + (event.cycle as f64 + phase_center) * period,
                    depth,
                    significance,
                    points: event.points,
                }
            })
            .collect();
        Some(Candidate {
            kind: CandidateKind::PeriodicDip,
            score,
            white_score: white,
            noise_factor,
            period,
            epoch: samples.origin + phase_center * period,
            duration: width as f64 / bins as f64 * period,
            depth,
            amplitude: depth,
            n_events,
            event_evidence,
            alias_ephemerides: Vec::new(),
            odd_even_sigma,
            note: if odd_even_sigma.is_some_and(|z| z >= 3.0) {
                "possible eclipsing binary: odd/even depth difference >= 3 sigma; vet blends and systematics".into()
            } else {
                "transit-like; vet blends, eclipses, and systematics".into()
            },
        })
    }
}

/// Insert `candidate` into the per-period shortlist: it replaces a weaker
/// window explaining the same events, and the list keeps the `limit` strongest.
fn keep_best(best: &mut Vec<Candidate>, candidate: Candidate, limit: usize) {
    if let Some(duplicate) = best
        .iter_mut()
        .find(|old| event_support_overlap(old, &candidate) >= 0.8)
    {
        if candidate.score > duplicate.score {
            *duplicate = candidate;
        }
    } else {
        best.push(candidate);
    }
    best.sort_by(|a, b| b.score.total_cmp(&a.score));
    best.truncate(limit);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{rednoise::RedNoise, segments::Segment};

    #[test]
    fn constant_residual_offset_does_not_change_transit_or_odd_even_fit() {
        let mut samples = Samples {
            origin: 0.0,
            times: (0..300).map(|i| i as f64 * 0.02).collect(),
            transit: (0..300)
                .map(|i| {
                    let phase = (i as f64 * 0.02).rem_euclid(1.0);
                    if (0.30..0.38).contains(&phase) {
                        -2000.0
                    } else {
                        0.0
                    }
                })
                .collect(),
            variability: vec![0.0; 300],
            transit_weights: vec![1.0 / 200.0_f64.powi(2); 300],
            variability_weights: vec![1.0; 300],
            cadence: 0.02,
            span: 5.98,
            segments: vec![Segment {
                start: 0,
                end: 300,
                cadence: 0.02,
            }],
            segment_ids: vec![0; 300],
            transit_noise: vec![200.0],
            red_noise: RedNoise::white(1),
        };
        let config = Config {
            min_period: 0.95,
            max_period: 1.05,
            periods: 100,
            durations: 6,
            threads: 1,
            ..Config::default()
        };
        let original = periodic_dips(&samples, config, &mut Vec::new()).candidates;
        for y in &mut samples.transit {
            *y += 5000.0;
        }
        let shifted = periodic_dips(&samples, config, &mut Vec::new()).candidates;
        assert!(!original.is_empty());
        assert_eq!(original.len(), shifted.len());
        assert!((original[0].score - shifted[0].score).abs() < 1e-9);
        assert_eq!(original[0].odd_even_sigma, shifted[0].odd_even_sigma);
    }

    /// Noise-free 4-day signal: a short deep dip and a broad shallow depression.
    fn two_dip_samples(red_noise: RedNoise) -> Samples {
        let cadence = 0.02;
        let count = 2400;
        let times: Vec<f64> = (0..count).map(|i| i as f64 * cadence).collect();
        let transit = times
            .iter()
            .map(|t| match (t / 4.0).rem_euclid(1.0) {
                p if (0.300..0.310).contains(&p) => -1500.0,
                p if (0.600..0.675).contains(&p) => -700.0,
                _ => 0.0,
            })
            .collect();
        Samples {
            origin: 0.0,
            span: times[count - 1],
            times,
            transit,
            variability: vec![0.0; count],
            transit_weights: vec![1.0 / 200.0_f64.powi(2); count],
            variability_weights: vec![1.0; count],
            cadence,
            segments: vec![Segment {
                start: 0,
                end: count,
                cadence,
            }],
            segment_ids: vec![0; count],
            transit_noise: vec![200.0],
            red_noise,
        }
    }

    fn config(scoring: Scoring) -> Config {
        Config {
            min_period: 3.9,
            max_period: 4.1,
            periods: 200,
            threads: 1,
            scoring,
            ..Config::default()
        }
    }

    #[test]
    fn beta_of_one_reproduces_the_white_statistic() {
        let samples = two_dip_samples(RedNoise::white(1));
        let white = periodic_dips(&samples, config(Scoring::White), &mut Vec::new()).candidates;
        let red = periodic_dips(&samples, config(Scoring::Correlated), &mut Vec::new()).candidates;
        assert!(!white.is_empty());
        assert_eq!(white.len(), red.len());
        for (w, r) in white.iter().zip(&red) {
            assert_eq!(
                (w.period, w.epoch, w.duration),
                (r.period, r.epoch, r.duration)
            );
            assert!((w.score / r.score - 1.0).abs() < 1e-12);
            assert_eq!(w.white_score, r.white_score);
            assert!(w.noise_factor.is_none());
            assert!((r.noise_factor.unwrap() - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn duration_dependent_beta_changes_which_window_is_selected() {
        // Beta^2 rises from 1 at 0.04 d to 9 at 0.35 d, so the broad 0.3-day
        // depression (higher white score) loses to the short dip.
        let red_noise = RedNoise {
            knots: vec![0.04, 0.35],
            beta2: vec![1.0, 9.0],
        };
        let samples = two_dip_samples(red_noise);
        let white = periodic_dips(&samples, config(Scoring::White), &mut Vec::new()).candidates;
        let red = periodic_dips(&samples, config(Scoring::Correlated), &mut Vec::new()).candidates;
        assert!(white[0].duration > 0.2, "{:?}", white[0]);
        assert!(red[0].duration < 0.1, "{:?}", red[0]);
        // The inflation is beta at the window's own duration (log-interpolated);
        // the out-of-window term keeps it slightly below.
        let beta = |d: f64| (1.0 + 8.0 * (d / 0.04).ln() / (0.35_f64 / 0.04).ln()).sqrt();
        let factor = red[0].noise_factor.unwrap();
        assert!(factor <= beta(red[0].duration) && factor > 0.99 * beta(red[0].duration));
        // The broad window is still reported, but its correlated score is lower.
        // A later search pass can pick a short slice of the depression, which
        // then carries the broad window as an alias sharing its events.
        let broad_score = red
            .iter()
            .find(|c| c.duration > 0.2)
            .map(|c| {
                assert!(c.noise_factor.unwrap() > 2.5, "{c:?}");
                c.score
            })
            .or_else(|| {
                red.iter()
                    .flat_map(|c| &c.alias_ephemerides)
                    .find(|a| a.duration > 0.2)
                    .map(|a| a.score)
            })
            .expect("broad window kept as a candidate or alias");
        assert!(broad_score < red[0].score);
        assert!(white[0].white_score > red[0].white_score);
    }

    #[test]
    fn periodogram_keeps_each_trial_best_and_candidate_limit_applies_after_ranking() {
        let samples = two_dip_samples(RedNoise::white(1));
        let all = Config {
            keep_periodogram: true,
            ..config(Scoring::Correlated)
        };
        let search = periodic_dips(&samples, all, &mut Vec::new());
        assert!(search.periodogram.len() > search.candidates.len());
        assert!(search
            .periodogram
            .windows(2)
            .all(|pair| pair[0].period < pair[1].period));
        let one = Config {
            transit_candidates: 1,
            ..config(Scoring::Correlated)
        };
        let limited = periodic_dips(&samples, one, &mut Vec::new()).candidates;
        assert_eq!(limited.len(), 1);
        assert_eq!(limited[0].score, search.candidates[0].score);
    }
}

//! Per-event (per-cycle) statistics of a phase window and the repeated-event vetting.

#[derive(Clone, Copy)]
pub(crate) struct CycleRange {
    pub cycle: i64,
    pub start: usize,
    pub end: usize,
}

#[derive(Clone, Copy)]
pub(crate) struct Epoch {
    pub cycle: i64,
    pub sum: f64,
    pub weight: f64,
    pub points: usize,
}

pub(crate) struct EpochIndex<'a> {
    pub cycles: &'a [CycleRange],
    pub point_bins: &'a [usize],
    pub point_sum_prefix: &'a [f64],
    pub point_weight_prefix: &'a [f64],
    /// Running count of positive-weight samples; masked samples do not count.
    pub point_count_prefix: &'a [usize],
}

/// Reconstruct event statistics only from phase-window points. Per-cycle
/// bins are sorted, so binary searches and prefix sums replace a full scan of
/// the light curve for every improving trial window.
impl EpochIndex<'_> {
    fn epoch_for_phase_range(
        &self,
        cycle: i64,
        range: CycleRange,
        lower_bin: usize,
        upper_bin: usize,
    ) -> Option<Epoch> {
        let bins = &self.point_bins[range.start..range.end];
        let start = range.start + bins.partition_point(|&bin| bin < lower_bin);
        let end = range.start + bins.partition_point(|&bin| bin < upper_bin);
        (end > start).then(|| Epoch {
            cycle,
            sum: self.point_sum_prefix[end] - self.point_sum_prefix[start],
            weight: self.point_weight_prefix[end] - self.point_weight_prefix[start],
            points: self.point_count_prefix[end] - self.point_count_prefix[start],
        })
    }

    pub(crate) fn collect_window_epochs(
        &self,
        start: usize,
        width: usize,
        bins: usize,
        epochs: &mut Vec<Epoch>,
        wrapped: &mut Vec<Epoch>,
        merged: &mut Vec<Epoch>,
    ) {
        epochs.clear();
        wrapped.clear();
        let end = start + width;
        if end <= bins {
            for &range in self.cycles {
                if let Some(epoch) = self.epoch_for_phase_range(range.cycle, range, start, end) {
                    epochs.push(epoch);
                }
            }
            return;
        }

        for &range in self.cycles {
            if let Some(epoch) = self.epoch_for_phase_range(range.cycle, range, start, bins) {
                epochs.push(epoch);
            }
            if let Some(epoch) = self.epoch_for_phase_range(range.cycle - 1, range, 0, end - bins) {
                wrapped.push(epoch);
            }
        }

        merged.clear();
        let (mut i, mut j) = (0, 0);
        while i < epochs.len() || j < wrapped.len() {
            match (epochs.get(i), wrapped.get(j)) {
                (Some(a), Some(b)) if a.cycle == b.cycle => {
                    merged.push(Epoch {
                        cycle: a.cycle,
                        sum: a.sum + b.sum,
                        weight: a.weight + b.weight,
                        points: a.points + b.points,
                    });
                    i += 1;
                    j += 1;
                }
                (Some(a), Some(b)) if a.cycle < b.cycle => {
                    merged.push(*a);
                    i += 1;
                }
                (Some(_), Some(b)) => {
                    merged.push(*b);
                    j += 1;
                }
                (Some(a), None) => {
                    merged.push(*a);
                    i += 1;
                }
                (None, Some(b)) => {
                    merged.push(*b);
                    j += 1;
                }
                (None, None) => break,
            }
        }
        std::mem::swap(epochs, merged);
    }
}

pub(crate) fn repeated_evidence(
    epochs: &[Epoch],
    out_mean: f64,
    out_weight: f64,
) -> Option<(usize, Option<f64>)> {
    let covered = epochs.iter().filter(|e| e.weight > 0.0).count();
    if covered < 2 {
        return None;
    }
    let positive_epochs = epochs
        .iter()
        .filter(|e| e.weight > 0.0 && e.points >= 2 && e.sum / e.weight < out_mean)
        .count();
    if positive_epochs < 2 || positive_epochs * 5 < covered * 3 {
        return None;
    }
    let contribution = |e: &Epoch| (out_mean * e.weight - e.sum).max(0.0);
    let largest = epochs
        .iter()
        .filter(|e| e.weight > 0.0)
        .map(contribution)
        .fold(0.0_f64, f64::max);
    if largest <= 0.0 {
        return None;
    }
    let (mut total, mut squares) = (0.0, 0.0);
    for epoch in epochs.iter().filter(|e| e.weight > 0.0) {
        let value = contribution(epoch) / largest;
        total += value;
        squares += value * value;
    }
    // Allow sub-3-sigma individual events to accumulate significance, but veto
    // a folded signal dominated by one isolated excursion. These thresholds
    // are screening heuristics, not calibrated detection probabilities.
    if total <= 1.25 || total * total / squares < (0.5 * covered as f64).max(1.5) {
        return None;
    }
    // A short fragment of one deep event can still inflate a broad box score.
    // Require corroboration after removing the strongest observed epoch.
    let strongest = epochs
        .iter()
        .enumerate()
        .filter(|(_, e)| e.weight > 0.0)
        .max_by(|(_, a), (_, b)| contribution(a).total_cmp(&contribution(b)))
        .unwrap()
        .0;
    let (mut remaining_sum, mut remaining_weight) = (0.0, 0.0);
    for (index, epoch) in epochs.iter().enumerate() {
        if index != strongest && epoch.weight > 0.0 {
            remaining_sum += epoch.sum;
            remaining_weight += epoch.weight;
        }
    }
    let corroboration = (out_mean - remaining_sum / remaining_weight)
        / (1.0 / remaining_weight + 1.0 / out_weight).sqrt();
    if !corroboration.is_finite() || corroboration < 3.0 {
        return None;
    }
    let mut parity = [(0.0_f64, 0.0_f64); 2];
    for epoch in epochs.iter().filter(|e| e.weight > 0.0) {
        let slot = epoch.cycle.rem_euclid(2) as usize;
        parity[slot].0 += epoch.sum;
        parity[slot].1 += epoch.weight;
    }
    let odd_even = if parity.iter().all(|v| v.1 > 0.0) {
        Some(
            (parity[0].0 / parity[0].1 - parity[1].0 / parity[1].1).abs()
                / (1.0 / parity[0].1 + 1.0 / parity[1].1).sqrt(),
        )
    } else {
        None
    };
    Some((positive_epochs, odd_even))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn one_deep_event_is_vetoed_but_weak_repeated_events_are_retained() {
        let repeated: Vec<Epoch> = (0..9)
            .map(|cycle| Epoch {
                cycle,
                sum: -2.0,
                weight: 1.0,
                points: 6,
            })
            .collect();
        assert_eq!(repeated_evidence(&repeated, 0.0, 1e12).unwrap().0, 9);
        let mut isolated = repeated;
        isolated[0].sum = -1000.0;
        assert!(repeated_evidence(&isolated, 0.0, 1e12).is_none());
    }

    #[test]
    fn unavailable_odd_even_statistic_is_not_reported_as_zero() {
        let epochs = vec![
            Epoch {
                cycle: 0,
                sum: -5.0,
                weight: 1.0,
                points: 3,
            },
            Epoch {
                cycle: 2,
                sum: -5.0,
                weight: 1.0,
                points: 3,
            },
        ];
        assert!(repeated_evidence(&epochs, 0.0, 1e12).unwrap().1.is_none());
    }

    #[test]
    fn indexed_epoch_collection_matches_direct_scan_for_wrapping_windows() {
        let point_cycles = [0_i64, 0, 1, 1, 2, 3, 3, 4];
        let point_bins = [2_usize, 9, 1, 7, 4, 0, 8, 3];
        let sums = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let weights = [1.0; 8];
        let mut sum_prefix = vec![0.0; sums.len() + 1];
        let mut weight_prefix = vec![0.0; sums.len() + 1];
        let mut cycles = Vec::<CycleRange>::new();
        for i in 0..sums.len() {
            sum_prefix[i + 1] = sum_prefix[i] + sums[i];
            weight_prefix[i + 1] = weight_prefix[i] + weights[i];
            if let Some(range) = cycles.last_mut().filter(|r| r.cycle == point_cycles[i]) {
                range.end = i + 1;
            } else {
                cycles.push(CycleRange {
                    cycle: point_cycles[i],
                    start: i,
                    end: i + 1,
                });
            }
        }
        let count_prefix: Vec<usize> = (0..=sums.len()).collect();
        let index = EpochIndex {
            cycles: &cycles,
            point_bins: &point_bins,
            point_sum_prefix: &sum_prefix,
            point_weight_prefix: &weight_prefix,
            point_count_prefix: &count_prefix,
        };
        for (start, width) in [(1, 4), (7, 4), (9, 2)] {
            let mut actual = Vec::new();
            index.collect_window_epochs(
                start,
                width,
                10,
                &mut actual,
                &mut Vec::new(),
                &mut Vec::new(),
            );
            let mut expected = Vec::<Epoch>::new();
            for i in 0..point_bins.len() {
                let b = point_bins[i];
                if (b + 10 - start) % 10 >= width {
                    continue;
                }
                let cycle = point_cycles[i] - i64::from(b < start);
                if let Some(last) = expected.last_mut().filter(|e| e.cycle == cycle) {
                    last.sum += sums[i];
                    last.weight += weights[i];
                    last.points += 1;
                } else {
                    expected.push(Epoch {
                        cycle,
                        sum: sums[i],
                        weight: weights[i],
                        points: 1,
                    });
                }
            }
            assert_eq!(actual.len(), expected.len());
            for (a, e) in actual.iter().zip(expected) {
                assert_eq!(a.cycle, e.cycle);
                assert_eq!(a.sum, e.sum);
                assert_eq!(a.weight, e.weight);
                assert_eq!(a.points, e.points);
            }
        }
    }
}

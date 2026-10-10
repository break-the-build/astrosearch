//! Deduplication and grouping of periodic-dip candidates, and alias labels.

use crate::{math::phase_distance, AliasEphemeris, AliasType, Candidate, CandidateKind};

pub(crate) fn collapse_periodic(
    mut candidates: Vec<Candidate>,
    span: f64,
    cap: usize,
) -> Vec<Candidate> {
    candidates.sort_by(|a, b| b.score.total_cmp(&a.score));
    let mut output: Vec<Candidate> = Vec::new();
    for candidate in candidates {
        let duplicate = output.iter().any(|old| {
            if candidate.kind != old.kind {
                return false;
            }
            if candidate.kind == CandidateKind::PeriodicDip {
                let drift =
                    (candidate.period - old.period).abs() / candidate.period.min(old.period) * span;
                drift <= candidate.duration.min(old.duration)
                    && phase_distance(candidate.epoch, old.epoch, old.period)
                        <= 0.5 * (candidate.duration + old.duration) + drift
            } else {
                (1.0 / candidate.period - 1.0 / old.period).abs() <= 0.5 / span
            }
        });
        if !duplicate {
            output.push(candidate);
        }
        if output.len() >= cap {
            break;
        }
    }
    output
}

/// Alternatives retained per event group.
const MAX_ALIASES: usize = 16;

/// Group transit ephemerides that explain the same observed events while
/// retaining their alternatives on the strongest representative.
pub(crate) fn collapse_transit_events(
    mut candidates: Vec<Candidate>,
    cap: usize,
) -> Vec<Candidate> {
    candidates.sort_by(|a, b| b.score.total_cmp(&a.score));
    let mut groups: Vec<Candidate> = Vec::new();
    for candidate in candidates {
        // Once every slot is taken and every group holds its maximum number
        // of alternatives, later candidates cannot change the output.
        if groups.len() >= cap
            && groups
                .iter()
                .all(|g| g.alias_ephemerides.len() >= MAX_ALIASES)
        {
            break;
        }
        let group = groups
            .iter_mut()
            .find(|primary| event_support_overlap(primary, &candidate) >= 0.8);
        if let Some(primary) = group {
            let distinct = (candidate.period - primary.period).abs() > f64::EPSILON
                || phase_distance(candidate.epoch, primary.epoch, primary.period)
                    > 0.5 * (candidate.duration + primary.duration);
            let event_span = candidate
                .event_evidence
                .first()
                .zip(candidate.event_evidence.last())
                .map_or(0.0, |(first, last)| last.epoch - first.epoch);
            let already_kept = primary.alias_ephemerides.iter().any(|alias| {
                (candidate.period - alias.period).abs() / candidate.period.min(alias.period)
                    * event_span
                    <= candidate.duration.min(alias.duration)
                    && phase_distance(candidate.epoch, alias.epoch, alias.period)
                        <= 0.5 * (candidate.duration + alias.duration)
            });
            if distinct && !already_kept && primary.alias_ephemerides.len() < MAX_ALIASES {
                primary.alias_ephemerides.push(AliasEphemeris {
                    period: candidate.period,
                    epoch: candidate.epoch,
                    duration: candidate.duration,
                    score: candidate.score,
                    n_events: candidate.n_events,
                    event_epochs: candidate
                        .event_evidence
                        .iter()
                        .map(|event| event.epoch)
                        .collect(),
                    alias_type: period_alias_type(primary.period, candidate.period),
                });
            }
        } else if groups.len() < cap {
            groups.push(candidate);
        }
    }
    groups
}

/// Harmonic relation of an alternative period to the representative's. The
/// longer period is compared with the shorter, within 1%: a ratio near 1 is
/// "same-period" (a shifted or resized window), P/2, P/3 and 2P are "integer",
/// and 2P/3 or 3P/2 are "half-integer".
fn period_alias_type(primary: f64, alias: f64) -> AliasType {
    let ratio = alias.max(primary) / alias.min(primary);
    let half_integer = (ratio * 2.0).round() / 2.0;
    if half_integer >= 1.0 && (ratio - half_integer).abs() / half_integer <= 0.01 {
        if half_integer == 1.0 {
            AliasType::SamePeriod
        } else if half_integer.fract() == 0.0 {
            AliasType::Integer
        } else {
            AliasType::HalfInteger
        }
    } else {
        AliasType::Unclassified
    }
}

/// Fraction of the smaller event set matched one-to-one, in time, by the
/// other. Event lists are in time order and separated by at least a period,
/// far more than the matching tolerance, so a two-pointer merge finds the
/// same matches as an exhaustive search in linear time.
pub(crate) fn event_support_overlap(a: &Candidate, b: &Candidate) -> f64 {
    if a.event_evidence.is_empty() || b.event_evidence.is_empty() {
        return 0.0;
    }
    let tolerance = 0.5 * (a.duration + b.duration);
    let (mut i, mut j, mut common) = (0, 0, 0);
    while i < a.event_evidence.len() && j < b.event_evidence.len() {
        let difference = a.event_evidence[i].epoch - b.event_evidence[j].epoch;
        if difference.abs() <= tolerance {
            common += 1;
            i += 1;
            j += 1;
        } else if difference < 0.0 {
            i += 1;
        } else {
            j += 1;
        }
    }
    common as f64 / a.event_evidence.len().min(b.event_evidence.len()) as f64
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::EventEvidence;

    fn candidate(period: f64, duration: f64, score: f64) -> Candidate {
        Candidate {
            kind: CandidateKind::PeriodicDip,
            score,
            white_score: score,
            noise_factor: None,
            period,
            epoch: 0.0,
            duration,
            depth: 1000.0,
            amplitude: 1000.0,
            n_events: 3,
            event_evidence: Vec::new(),
            alias_ephemerides: Vec::new(),
            odd_even_sigma: None,
            note: String::new(),
        }
    }

    fn with_events(mut candidate: Candidate, epochs: &[f64]) -> Candidate {
        candidate.event_evidence = epochs
            .iter()
            .map(|&epoch| EventEvidence {
                epoch,
                depth: 500.0,
                significance: 5.0,
                points: 4,
            })
            .collect();
        candidate
    }

    /// Exhaustive one-to-one matching, the previous implementation.
    fn exhaustive_overlap(a: &Candidate, b: &Candidate) -> f64 {
        let tolerance = 0.5 * (a.duration + b.duration);
        let mut matched = vec![false; b.event_evidence.len()];
        let mut common = 0;
        for event in &a.event_evidence {
            if let Some(index) = (0..b.event_evidence.len()).find(|&k| {
                !matched[k] && (event.epoch - b.event_evidence[k].epoch).abs() <= tolerance
            }) {
                matched[index] = true;
                common += 1;
            }
        }
        common as f64 / a.event_evidence.len().min(b.event_evidence.len()) as f64
    }

    #[test]
    fn event_grouping_keeps_alias_ephemerides_without_spending_candidate_slots() {
        let output = collapse_transit_events(
            vec![
                with_events(candidate(3.0, 0.1, 20.0), &[1.0, 11.0]),
                with_events(candidate(3.001, 0.1, 19.0), &[1.02, 11.01]),
                with_events(candidate(4.0, 0.1, 18.0), &[2.0, 14.0]),
            ],
            2,
        );
        assert_eq!(output.len(), 2);
        assert_eq!(output[0].period, 3.0);
        assert_eq!(output[0].alias_ephemerides.len(), 1);
        assert!((output[0].alias_ephemerides[0].period - 3.001).abs() < 1e-12);
        assert_eq!(output[1].period, 4.0);
    }

    #[test]
    fn event_grouping_preserves_distinct_planets_that_share_one_event() {
        let output = collapse_transit_events(
            vec![
                with_events(candidate(3.0, 0.1, 20.0), &[1.0, 11.0]),
                with_events(candidate(4.0, 0.1, 19.0), &[1.02, 15.0]),
            ],
            12,
        );
        assert_eq!(output.len(), 2);
        assert!(output
            .iter()
            .all(|candidate| candidate.alias_ephemerides.is_empty()));
    }

    #[test]
    fn event_grouping_labels_harmonic_solutions_as_aliases() {
        let output = collapse_transit_events(
            vec![
                with_events(candidate(6.0, 0.1, 20.0), &[1.0, 7.0]),
                with_events(candidate(3.0, 0.1, 19.0), &[1.0, 4.0, 7.0]),
            ],
            12,
        );
        assert_eq!(output.len(), 1);
        assert_eq!(output[0].alias_ephemerides.len(), 1);
        assert_eq!(output[0].alias_ephemerides[0].period, 3.0);
        assert_eq!(output[0].alias_ephemerides[0].alias_type, "integer");
    }

    #[test]
    fn alias_labels_cover_harmonics_and_subharmonics() {
        for (alias, label) in [
            (3.02, "same-period"),
            (2.98, "same-period"),
            (6.0, "integer"),
            (1.5, "integer"),
            (1.0, "integer"),
            (4.5, "half-integer"),
            (2.0, "half-integer"),
            (3.7, "unclassified"),
        ] {
            assert_eq!(period_alias_type(3.0, alias), label, "alias {alias}");
        }
    }

    #[test]
    fn linear_overlap_matches_exhaustive_matching() {
        let mut state = 9_u64;
        let mut uniform = || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (state >> 11) as f64 / (1_u64 << 53) as f64
        };
        for _ in 0..500 {
            let events = |period: f64, epoch: f64, count: usize, keep: &mut dyn FnMut() -> f64| {
                (0..count)
                    .map(|k| epoch + k as f64 * period)
                    .filter(|_| keep() > 0.3)
                    .collect::<Vec<f64>>()
            };
            let (pa, pb) = (1.0 + 9.0 * uniform(), 1.0 + 9.0 * uniform());
            let (ea, eb) = (pa * uniform(), pb * uniform());
            let a = with_events(candidate(pa, 0.1, 10.0), &events(pa, ea, 30, &mut uniform));
            let b = with_events(candidate(pb, 0.15, 9.0), &events(pb, eb, 30, &mut uniform));
            if a.event_evidence.is_empty() || b.event_evidence.is_empty() {
                continue;
            }
            assert_eq!(event_support_overlap(&a, &b), exhaustive_overlap(&a, &b));
            assert_eq!(event_support_overlap(&b, &a), exhaustive_overlap(&b, &a));
        }
    }

    #[test]
    fn grouping_stops_early_only_when_output_cannot_change() {
        let mut input = vec![with_events(candidate(3.0, 0.1, 100.0), &[1.0, 4.0, 7.0])];
        // Twenty distinct aliases of the same events (each drifts more than a
        // duration from the others over the 6-day event span); 16 are kept.
        for k in 0..20 {
            let period = 3.0 + 0.2 * (k + 1) as f64;
            input.push(with_events(
                candidate(period, 0.1, 99.0 - k as f64),
                &[1.0, 4.0, 7.0],
            ));
        }
        input.push(with_events(candidate(5.0, 0.1, 1.0), &[2.0, 12.0]));
        let one = collapse_transit_events(input.clone(), 1);
        assert_eq!(one.len(), 1);
        assert_eq!(one[0].alias_ephemerides.len(), MAX_ALIASES);
        let two = collapse_transit_events(input, 2);
        assert_eq!(two.len(), 2);
        assert_eq!(two[1].period, 5.0);
    }

    #[test]
    fn nearby_resolved_periods_survive_and_near_identical_trials_merge() {
        let output = collapse_periodic(
            vec![
                candidate(3.0, 0.1, 20.0),
                candidate(3.00001, 0.1, 19.0),
                candidate(3.2, 0.1, 18.0),
            ],
            27.0,
            12,
        );
        assert_eq!(output.len(), 2);
        assert_eq!(output[1].period, 3.2);
    }
}

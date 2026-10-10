//! Split a light curve into stretches with one cadence and a shared noise level,
//! such as a TESS sector (or orbit) observed with one data product.

use crate::math::median;

/// Cadences within this factor are treated as the same observing mode.
pub(crate) const CADENCE_RATIO: f64 = 1.5;
/// Samples that establish, and confirm a change of, the native cadence. A switch
/// must persist this long, so quality-masked stretches stay in their sector.
const SWITCH_WINDOW: usize = 50;
/// Quality masking removes samples and widens some spacings, but any surviving
/// consecutive pair still shows the native cadence; a low quantile tracks it.
const NATIVE_QUANTILE: f64 = 0.2;
/// Same-cadence pieces closer than this many cadences stay one segment.
const GAP_CADENCES: f64 = 5.0;

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct Segment {
    pub start: usize,
    /// Exclusive.
    pub end: usize,
    /// Median spacing inside the segment; the global cadence for one-point segments.
    pub cadence: f64,
}

impl Segment {
    pub fn len(&self) -> usize {
        self.end - self.start
    }
}

pub(crate) fn ratio(a: f64, b: f64) -> f64 {
    a.max(b) / a.min(b)
}

/// Boundaries fall at gaps longer than both `width` and five segment cadences,
/// and where the sampling cadence changes by more than `CADENCE_RATIO`.
/// `times` must be sorted and distinct.
pub(crate) fn split(times: &[f64], width: f64, fallback_cadence: f64) -> Vec<Segment> {
    let mut pieces = Vec::new();
    let mut start = 0;
    while start < times.len() {
        let mut end = start + 1;
        while end < times.len() && times[end] - times[end - 1] <= width {
            end += 1;
        }
        split_cadence_changes(times, start, end, &mut pieces);
        start = end;
    }
    let pieces: Vec<Segment> = pieces
        .into_iter()
        .map(|(start, end)| with_cadence(times, start, end, fallback_cadence))
        .collect();
    merge_same_cadence(times, pieces, width, fallback_cadence)
}

/// Splits one gap-free run where a different spacing persists.
fn split_cadence_changes(
    times: &[f64],
    start: usize,
    end: usize,
    pieces: &mut Vec<(usize, usize)>,
) {
    // Each point's spacing is its smaller neighbor difference within the run,
    // so the first point after a switch already reports the new cadence.
    let spacing: Vec<f64> = (start..end)
        .map(|i| {
            let before = (i > start).then(|| times[i] - times[i - 1]);
            let after = (i + 1 < end).then(|| times[i + 1] - times[i]);
            match (before, after) {
                (Some(a), Some(b)) => a.min(b),
                (Some(a), None) | (None, Some(a)) => a,
                (None, None) => f64::NAN,
            }
        })
        .collect();
    let window = |from: usize| native(&spacing[from..(from + SWITCH_WINDOW).min(spacing.len())]);
    let mut piece_start = 0;
    let mut reference = window(0);
    // A change in the last SWITCH_WINDOW samples cannot be confirmed and stays.
    let last = spacing.len().saturating_sub(SWITCH_WINDOW - 1);
    for (j, &value) in spacing.iter().enumerate().take(last).skip(1) {
        if ratio(value, reference) <= CADENCE_RATIO {
            continue;
        }
        let confirmed = window(j);
        if ratio(confirmed, reference) > CADENCE_RATIO {
            pieces.push((start + piece_start, start + j));
            piece_start = j;
            reference = confirmed;
        }
    }
    pieces.push((start + piece_start, end));
}

/// Low-quantile spacing; NaN when no finite spacing is available.
fn native(spacing: &[f64]) -> f64 {
    let mut values: Vec<f64> = spacing.iter().copied().filter(|v| v.is_finite()).collect();
    if values.is_empty() {
        return f64::NAN;
    }
    let k = ((values.len() - 1) as f64 * NATIVE_QUANTILE) as usize;
    *values.select_nth_unstable_by(k, f64::total_cmp).1
}

fn with_cadence(times: &[f64], start: usize, end: usize, fallback: f64) -> Segment {
    let cadence = if end - start > 1 {
        median(times[start..end].windows(2).map(|t| t[1] - t[0]))
    } else {
        fallback
    };
    Segment {
        start,
        end,
        cadence,
    }
}

/// Coarse cadences can exceed `width`; keep their ordinary spacing together.
fn merge_same_cadence(
    times: &[f64],
    pieces: Vec<Segment>,
    width: f64,
    fallback: f64,
) -> Vec<Segment> {
    let mut merged: Vec<Segment> = Vec::with_capacity(pieces.len());
    for piece in pieces {
        let joins = merged.last().is_some_and(|last| {
            let gap = times[piece.start] - times[last.end - 1];
            ratio(last.cadence, piece.cadence) <= CADENCE_RATIO
                && gap <= width.max(GAP_CADENCES * last.cadence.max(piece.cadence))
        });
        if joins {
            let last = merged.pop().unwrap();
            merged.push(with_cadence(times, last.start, piece.end, fallback));
        } else {
            merged.push(piece);
        }
    }
    merged
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run(start: f64, cadence: f64, count: usize) -> Vec<f64> {
        (0..count).map(|i| start + i as f64 * cadence).collect()
    }

    fn bounds(segments: &[Segment]) -> Vec<(usize, usize)> {
        segments.iter().map(|s| (s.start, s.end)).collect()
    }

    #[test]
    fn sectors_separated_by_long_gaps_keep_their_own_cadence() {
        let mut times = run(0.0, 1800.0 / 86400.0, 100);
        times.extend(run(30.0, 120.0 / 86400.0, 500));
        let segments = split(&times, 0.5, 120.0 / 86400.0);
        assert_eq!(bounds(&segments), vec![(0, 100), (100, 600)]);
        assert!((segments[0].cadence * 86400.0 - 1800.0).abs() < 1e-6);
        assert!((segments[1].cadence * 86400.0 - 120.0).abs() < 1e-6);
    }

    #[test]
    fn cadence_switch_without_a_gap_splits_in_both_directions() {
        let fast = 120.0 / 86400.0;
        let slow = 1800.0 / 86400.0;
        let mut times = run(0.0, slow, 80);
        times.extend(run(times[79] + fast, fast, 200));
        times.extend(run(times[279] + slow, slow, 80));
        let segments = split(&times, 0.5, fast);
        assert_eq!(segments.len(), 3, "{segments:?}");
        // The last slow sample's smaller spacing is the following fast one.
        assert!((78..=80).contains(&segments[0].end));
        assert!((279..=281).contains(&segments[1].end));
        assert!(ratio(segments[0].cadence, slow) < 1.01);
        assert!(ratio(segments[1].cadence, fast) < 1.01);
        assert!(ratio(segments[2].cadence, slow) < 1.01);
    }

    #[test]
    fn short_dropouts_and_isolated_samples_do_not_fragment_a_sector() {
        let cadence = 120.0 / 86400.0;
        let mut times = run(0.0, cadence, 300);
        times.retain(|t| !(0.1..0.15).contains(t) && !(0.2..0.21).contains(t));
        let segments = split(&times, 0.5, cadence);
        assert_eq!(bounds(&segments), vec![(0, times.len())]);
    }

    #[test]
    fn quality_masked_stretches_stay_in_their_sector() {
        // Modeled on QLP sector 98: inside a 200-second sector, stretches lose
        // every other sample, or 2 of every 3, leaving 400-1000 s spacings.
        let cadence = 200.0 / 86400.0;
        let times: Vec<f64> = run(0.0, cadence, 3000)
            .into_iter()
            .enumerate()
            .filter(|(i, _)| match i {
                1000..1030 => i % 2 == 0,
                1500..1800 => i % 7 < 3,
                2200..2215 => i % 5 == 0,
                _ => true,
            })
            .map(|(_, t)| t)
            .collect();
        let segments = split(&times, 0.5, cadence);
        assert_eq!(bounds(&segments), vec![(0, times.len())], "{segments:?}");
        assert!(ratio(segments[0].cadence, cadence) < 1.01);
    }

    #[test]
    fn orbit_gap_splits_and_coarse_cadence_matches_previous_gap_rule() {
        let cadence = 120.0 / 86400.0;
        let mut times = run(0.0, cadence, 200);
        times.extend(run(1.5, cadence, 200));
        assert_eq!(split(&times, 0.5, cadence).len(), 2);
        // A 0.2-day cadence keeps a one-sample (0.4-day) gap together, as
        // max(width, 5 cadences) did before segmentation was introduced.
        let mut coarse = run(0.0, 0.2, 20);
        coarse.extend(run(4.2, 0.2, 20));
        assert_eq!(split(&coarse, 0.3, 0.2).len(), 1);
    }
}

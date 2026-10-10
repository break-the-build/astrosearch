//! Input samples and the analysis results returned to callers.

use std::fmt;

/// Implements `as_str`, `Display`, and comparison with string literals, so
/// the CSV text and existing `kind == "periodic_dip"` comparisons are unchanged.
macro_rules! labeled_enum {
    ($name:ident { $($variant:ident => $label:literal),+ $(,)? }) => {
        impl $name {
            pub fn as_str(self) -> &'static str {
                match self {
                    $(Self::$variant => $label),+
                }
            }
        }
        impl fmt::Display for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                f.write_str(self.as_str())
            }
        }
        impl PartialEq<&str> for $name {
            fn eq(&self, other: &&str) -> bool {
                self.as_str() == *other
            }
        }
    };
}

/// Detector that produced a candidate.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CandidateKind {
    /// Box-shaped dips repeating in three or more observed transits.
    PeriodicDip,
    /// Exactly two observed transits: the period is ambiguous (any integer
    /// fraction of their separation that predicts no contradicting data), so
    /// it needs alias enumeration rather than a single ephemeris.
    DuoTransit,
    PeriodicVariability,
    IsolatedDip,
    FlareOrBrightening,
}

labeled_enum!(CandidateKind {
    PeriodicDip => "periodic_dip",
    DuoTransit => "duo_transit",
    PeriodicVariability => "periodic_variability",
    IsolatedDip => "isolated_dip",
    FlareOrBrightening => "flare_or_brightening",
});

impl CandidateKind {
    /// Transit-like candidates from the box search (periodic or duo).
    pub fn is_transit(self) -> bool {
        matches!(self, Self::PeriodicDip | Self::DuoTransit)
    }
}

/// Harmonic relation of an alternative ephemeris to its group representative:
/// the longer of the two periods divided by the shorter, within 1%.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AliasType {
    /// Ratio near 1: a shifted or resized window at essentially the same period.
    SamePeriod,
    /// Ratio near 2, 3, ... (for example 2P, P/2, P/3).
    Integer,
    /// Ratio near 1.5, 2.5, ... (for example 3P/2, 2P/3).
    HalfInteger,
    Unclassified,
}

labeled_enum!(AliasType {
    SamePeriod => "same-period",
    Integer => "integer",
    HalfInteger => "half-integer",
    Unclassified => "unclassified",
});

#[derive(Clone, Copy, Debug)]
pub struct Point {
    pub time: f64,
    pub flux: f64,
    /// Uncertainty in the same units as flux; zero means unavailable.
    pub error: f64,
}

#[derive(Clone, Debug)]
pub struct Candidate {
    pub kind: CandidateKind,
    /// The ranking statistic (see `Scoring`), with no survey-wide false-alarm calibration.
    pub score: f64,
    /// The original independent-error statistic of the same window.
    pub white_score: f64,
    /// `white_score / score` for periodic dips searched with correlated scoring:
    /// the window's effective correlated-noise inflation. None otherwise.
    pub noise_factor: Option<f64>,
    pub period: f64,
    /// Center of the fitted transit window or observed isolated excursion.
    pub epoch: f64,
    pub duration: f64,
    pub depth: f64,
    pub amplitude: f64,
    pub n_events: usize,
    /// Per-transit observed support for periodic-dip candidates; empty otherwise.
    pub event_evidence: Vec<EventEvidence>,
    /// Alternative ephemerides grouped with this candidate in event-grouping mode.
    pub alias_ephemerides: Vec<AliasEphemeris>,
    pub odd_even_sigma: Option<f64>,
    pub note: String,
}

impl Candidate {
    /// Fraction of the event-level signal (sum of positive event S/N squared)
    /// carried by the single strongest event; 1 without event evidence.
    pub fn single_event_share(&self) -> f64 {
        let squares = self
            .event_evidence
            .iter()
            .map(|e| e.significance.max(0.0).powi(2));
        let (total, largest) = squares.fold((0.0, 0.0_f64), |(t, m), s| (t + s, m.max(s)));
        if total > 0.0 {
            largest / total
        } else {
            1.0
        }
    }

    /// Score left after removing the strongest event, `score * sqrt(1 - share)`
    /// (score^2 approximately sums event S/N^2). On the known-answer benchmark,
    /// null-control maxima are dominated by one event (median share 0.60 vs
    /// 0.11 for planets), so this separates planets from systematics better
    /// than `score` (calibration-split AUC 0.91 vs 0.80).
    pub fn jackknife_score(&self) -> f64 {
        self.score * (1.0 - self.single_event_share()).max(0.0).sqrt()
    }
}

#[derive(Clone, Debug)]
pub struct EventEvidence {
    /// Predicted midpoint of this observed event in the input time system.
    pub epoch: f64,
    /// Measured local dip relative to the candidate's out-of-transit baseline, ppm.
    pub depth: f64,
    /// Independent-error event S/N diagnostic; not a calibrated probability.
    pub significance: f64,
    pub points: usize,
}

#[derive(Clone, Debug)]
pub struct AliasEphemeris {
    pub period: f64,
    pub epoch: f64,
    pub duration: f64,
    pub score: f64,
    pub n_events: usize,
    /// Midpoints of observed events supporting this alternative solution.
    pub event_epochs: Vec<f64>,
    pub alias_type: AliasType,
}

/// One independent stretch with a single cadence, such as a TESS sector or orbit.
#[derive(Clone, Debug)]
pub struct SegmentSummary {
    pub start_time: f64,
    pub end_time: f64,
    pub observations: usize,
    pub cadence_days: f64,
    /// Robust neighbor-difference scatter of the detrended flux.
    pub noise_ppm: f64,
    /// Correlated-noise factor beta at each duration knot, as (days, beta).
    pub red_noise: Vec<(f64, f64)>,
}

#[derive(Debug)]
pub struct Report {
    pub observations: usize,
    pub span_days: f64,
    /// Number of trial periods evaluated by the transit detector.
    pub transit_trials: usize,
    pub segments: Vec<SegmentSummary>,
    pub candidates: Vec<Candidate>,
    /// Best accepted window per trial period, in period order, when requested.
    pub periodogram: Vec<Candidate>,
    pub warnings: Vec<String>,
}

#[cfg(test)]
mod tests {
    use super::*;

    fn candidate(significances: &[f64]) -> Candidate {
        Candidate {
            kind: CandidateKind::PeriodicDip,
            score: 30.0,
            white_score: 30.0,
            noise_factor: None,
            period: 3.0,
            epoch: 1.0,
            duration: 0.1,
            depth: 1000.0,
            amplitude: 0.0,
            n_events: significances.len(),
            event_evidence: significances
                .iter()
                .map(|&significance| EventEvidence {
                    epoch: 1.0,
                    depth: 1000.0,
                    significance,
                    points: 10,
                })
                .collect(),
            alias_ephemerides: Vec::new(),
            odd_even_sigma: None,
            note: String::new(),
        }
    }

    #[test]
    fn jackknife_score_discounts_single_event_dominance() {
        let even = candidate(&[5.0; 9]);
        assert!((even.single_event_share() - 1.0 / 9.0).abs() < 1e-12);
        assert!((even.jackknife_score() - 30.0 * (8.0_f64 / 9.0).sqrt()).abs() < 1e-9);
        let one = candidate(&[20.0, 1.0, -2.0, 1.0]);
        assert!(one.single_event_share() > 0.99);
        assert!(one.jackknife_score() < 3.0);
        assert_eq!(candidate(&[]).jackknife_score(), 0.0);
    }
}

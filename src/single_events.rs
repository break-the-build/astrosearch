//! Isolated (non-periodic) brightening and dimming events.

use crate::{
    math::{median, phase_distance},
    preprocess::Samples,
    Candidate, CandidateKind,
};

pub(crate) fn events(samples: &Samples, min_points: usize) -> Vec<Candidate> {
    let baseline = median(samples.transit.iter().copied());
    let z: Vec<f64> = samples
        .transit
        .iter()
        .zip(&samples.transit_weights)
        .map(|(y, w)| (y - baseline) * w.sqrt())
        .collect();
    let mut candidates = Vec::new();
    let mut i = 0;
    while i < z.len() {
        // Include both shoulders of an excursion; require a 4.5-sigma peak
        // somewhere in the contiguous, same-sign 2.5-sigma run.
        if z[i].abs() < 2.5 {
            i += 1;
            continue;
        }
        let sign = z[i].signum();
        let start = i;
        // Runs stay inside one cadence segment and break at gaps of its own cadence.
        let segment = samples.segment_ids[start];
        let cadence = samples.segments[segment].cadence;
        let (mut weighted_sum, mut weight, mut peak_z, mut peak) = (0.0, 0.0, 0.0_f64, 0.0_f64);
        while i < z.len()
            && z[i] * sign >= 2.5
            && (i == start
                || (samples.segment_ids[i] == segment
                    && samples.times[i] - samples.times[i - 1] <= 3.0 * cadence))
        {
            let residual = samples.transit[i] - baseline;
            weighted_sum += residual * samples.transit_weights[i];
            weight += samples.transit_weights[i];
            peak_z = peak_z.max(z[i] * sign);
            peak = peak.max(residual * sign);
            i += 1;
        }
        if i - start < min_points || peak_z < 4.5 || weight <= 0.0 {
            continue;
        }
        let score = weighted_sum.abs() / weight.sqrt();
        candidates.push(Candidate {
            kind: if sign > 0.0 {
                CandidateKind::FlareOrBrightening
            } else {
                CandidateKind::IsolatedDip
            },
            score,
            white_score: score,
            noise_factor: None,
            period: 0.0,
            epoch: samples.origin + 0.5 * (samples.times[start] + samples.times[i - 1]),
            duration: samples.times[i - 1] - samples.times[start] + cadence,
            depth: if sign < 0.0 {
                -weighted_sum / weight
            } else {
                0.0
            },
            amplitude: peak,
            n_events: 1,
            event_evidence: Vec::new(),
            alias_ephemerides: Vec::new(),
            odd_even_sigma: None,
            note: "single-event trigger; inspect pixels and nearby sources".into(),
        });
    }
    candidates.sort_by(|a, b| b.score.total_cmp(&a.score));
    // Runs are already disjoint; nearby opposite-sign events remain independent.
    candidates.truncate(100);
    candidates
}

pub(crate) fn suppress_periodic_event_duplicates(candidates: &mut Vec<Candidate>, cadence: f64) {
    let periodic: Vec<(f64, f64, f64)> = candidates
        .iter()
        .filter(|c| c.kind.is_transit() && c.score >= 8.0)
        .map(|c| (c.epoch, c.period, c.duration))
        .collect();
    candidates.retain(|event| {
        event.kind != CandidateKind::IsolatedDip
            || !periodic.iter().any(|(epoch, period, duration)| {
                // Both epochs describe window centers. Suppress only contained events,
                // rather than dropping a separate excursion merely near a transit.
                phase_distance(event.epoch, *epoch, *period) + event.duration / 2.0
                    <= duration / 2.0 + cadence
            })
    });
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{rednoise::RedNoise, segments::Segment};

    #[test]
    fn event_shoulders_and_uncertainties_define_complete_window() {
        let mut samples = Samples {
            origin: 100.0,
            times: (0..30).map(|i| i as f64).collect(),
            transit: vec![0.0; 30],
            variability: vec![0.0; 30],
            transit_weights: vec![1.0; 30],
            variability_weights: vec![1.0; 30],
            cadence: 1.0,
            span: 29.0,
            segments: vec![Segment {
                start: 0,
                end: 30,
                cadence: 1.0,
            }],
            segment_ids: vec![0; 30],
            transit_noise: vec![1.0],
            red_noise: RedNoise::white(1),
        };
        samples.transit[10..13].copy_from_slice(&[3.0, 6.0, 3.0]);
        let output = events(&samples, 2);
        assert_eq!(output.len(), 1);
        assert_eq!(output[0].epoch, 111.0);
        assert_eq!(output[0].duration, 3.0);
        samples.transit_weights[10..13].fill(0.01);
        assert!(events(&samples, 2).is_empty());
    }

    #[test]
    fn event_runs_use_their_segment_cadence_and_stop_at_segment_boundaries() {
        // 40 fast samples set the global median cadence; 10 slow samples follow.
        let mut times: Vec<f64> = (0..40).map(|i| i as f64).collect();
        times.extend((0..10).map(|i| 40.0 + i as f64 * 15.0));
        let mut samples = Samples {
            origin: 0.0,
            span: *times.last().unwrap(),
            times,
            transit: vec![0.0; 50],
            variability: vec![0.0; 50],
            transit_weights: vec![1.0; 50],
            variability_weights: vec![1.0; 50],
            cadence: 1.0,
            segments: vec![
                Segment {
                    start: 0,
                    end: 40,
                    cadence: 1.0,
                },
                Segment {
                    start: 40,
                    end: 50,
                    cadence: 15.0,
                },
            ],
            segment_ids: [vec![0; 40], vec![1; 10]].concat(),
            transit_noise: vec![1.0, 1.0],
            red_noise: RedNoise::white(2),
        };
        // Two consecutive slow samples are 15 global cadences apart.
        samples.transit[44..46].copy_from_slice(&[-5.0, -5.0]);
        let output = events(&samples, 2);
        assert_eq!(output.len(), 1, "{output:?}");
        assert_eq!(output[0].kind, "isolated_dip");
        assert_eq!(output[0].duration, 30.0);
        // An excursion spanning the boundary is two runs, each too short.
        samples.transit[44..46].fill(0.0);
        samples.transit[39..41].copy_from_slice(&[-5.0, -5.0]);
        assert!(events(&samples, 2).is_empty());
    }
}

//! Correlated-noise scoring on mixed-cadence light curves.

use astrosearch::{
    analyze, render_csv, Candidate, CandidateGrouping, Config, Point, Report, Scoring,
};

struct Noise(u64);

impl Noise {
    fn uniform(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (((self.0 >> 11) as f64) + 0.5) / ((1_u64 << 53) as f64)
    }
    fn gaussian(&mut self) -> f64 {
        (-2.0 * self.uniform().ln()).sqrt() * (std::f64::consts::TAU * self.uniform()).cos()
    }
}

/// A white 30-minute sector, then a 10-minute sector 60 days later with an
/// Ornstein-Uhlenbeck process (0.12 d) and two dip-like glitches per orbit that
/// the quoted errors omit. Each sector has two 10-day orbits.
fn curve(seed: u64, red: bool) -> Vec<Point> {
    let mut noise = Noise(seed);
    let mut points = Vec::new();
    for (start, seconds, white, ou) in [(0.0, 1800.0, 250e-6, 0.0), (60.0, 600.0, 433e-6, 170e-6)] {
        let step = seconds / 86400.0;
        let decay = (-step / 0.12_f64).exp();
        for orbit in [start, start + 11.0] {
            let glitches: Vec<f64> = (0..2).map(|_| orbit + 10.0 * noise.uniform()).collect();
            let mut state = 0.0;
            for i in 0..(10.0 / step) as usize {
                let time = orbit + i as f64 * step;
                let mut flux = 1.0 + white * noise.gaussian();
                if red && ou > 0.0 {
                    state = state * decay + ou * (1.0 - decay * decay).sqrt() * noise.gaussian();
                    flux += state;
                    for g in &glitches {
                        flux -= 250e-6 * (-0.5 * ((time - g) / 0.04).powi(2)).exp();
                    }
                }
                points.push(Point {
                    time,
                    flux,
                    error: white,
                });
            }
        }
    }
    points
}

fn config(scoring: Scoring, threads: usize) -> Config {
    Config {
        min_period: 1.0,
        max_period: 10.0,
        scoring,
        threads,
        ..Config::default()
    }
}

fn dips(report: &Report) -> Vec<&Candidate> {
    report
        .candidates
        .iter()
        .filter(|c| c.kind.is_transit())
        .collect()
}

#[test]
fn red_segment_reports_larger_beta_than_white_segment() {
    let report = analyze(curve(1, true), config(Scoring::Correlated, 0)).unwrap();
    assert_eq!(report.segments.len(), 4);
    let beta_at = |segment: usize, days: f64| {
        report.segments[segment]
            .red_noise
            .iter()
            .min_by(|a, b| (a.0 / days).ln().abs().total_cmp(&(b.0 / days).ln().abs()))
            .unwrap()
            .1
    };
    for white in [0, 1] {
        assert!(
            beta_at(white, 2.0 / 24.0) < 1.12,
            "{:?}",
            report.segments[white]
        );
    }
    for red in [2, 3] {
        assert!(beta_at(red, 2.0 / 24.0) > 1.2, "{:?}", report.segments[red]);
        assert!(beta_at(red, 2.0 / 24.0) > beta_at(red, 0.01));
    }
}

#[test]
fn correlated_scoring_suppresses_red_noise_false_alarms_at_the_default_threshold() {
    let (mut white, mut correlated) = (0, 0);
    for seed in 0..8 {
        let points = curve(100 + seed, true);
        let w = analyze(points.clone(), config(Scoring::White, 0)).unwrap();
        let r = analyze(points, config(Scoring::Correlated, 0)).unwrap();
        white += usize::from(!dips(&w).is_empty());
        correlated += usize::from(!dips(&r).is_empty());
        // Every candidate keeps both statistics; correlated never exceeds white.
        for c in dips(&r) {
            assert!(
                c.score <= c.white_score && c.noise_factor.unwrap() >= 1.0,
                "{c:?}"
            );
        }
    }
    // This small deterministic fixture checks relative suppression, not the
    // absolute false-alarm rate. Eight realizations are insufficient to require
    // a particular baseline count; calibration belongs in the larger harness.
    assert!(
        white > 0,
        "fixture produced no white-score triggers to suppress"
    );
    assert!(
        correlated * 2 <= white,
        "correlated {correlated} vs white {white}"
    );
}

#[test]
fn injected_transit_in_red_noise_is_still_recovered() {
    let mut points = curve(7, true);
    let (period, epoch, depth) = (3.37, 1.2, 600e-6);
    for p in &mut points {
        let phase = (p.time - epoch + period / 2.0).rem_euclid(period) - period / 2.0;
        if phase.abs() < 1.0 / 24.0 {
            p.flux -= depth;
        }
    }
    let report = analyze(points, config(Scoring::Correlated, 0)).unwrap();
    let best = dips(&report)[0];
    assert!(
        (best.period - period).abs() * report.span_days / period < 2.0 / 24.0,
        "{best:?}"
    );
    assert!(best.noise_factor.unwrap() > 1.0);
}

#[test]
fn correlated_scoring_is_bit_identical_serial_and_parallel() {
    let mut points = curve(21, true);
    for p in &mut points {
        if ((p.time - 0.7).rem_euclid(2.9)) < 0.08 {
            p.flux -= 500e-6;
        }
    }
    let base = Config {
        keep_periodogram: true,
        ..config(Scoring::Correlated, 1)
    };
    let serial = analyze(points.clone(), base).unwrap();
    assert!(!serial.periodogram.is_empty());
    let bits = |c: &Candidate| {
        [
            c.score,
            c.white_score,
            c.period,
            c.epoch,
            c.duration,
            c.depth,
        ]
        .map(f64::to_bits)
        .to_vec()
        .into_iter()
        .chain(c.noise_factor.map(f64::to_bits))
        .collect::<Vec<_>>()
    };
    for threads in [4, 24] {
        let parallel = analyze(points.clone(), Config { threads, ..base }).unwrap();
        for (a, b) in [
            (&serial.candidates, &parallel.candidates),
            (&serial.periodogram, &parallel.periodogram),
        ] {
            assert_eq!(a.len(), b.len());
            assert!(a.iter().zip(b).all(|(x, y)| bits(x) == bits(y)));
        }
        for (a, b) in serial.segments.iter().zip(&parallel.segments) {
            assert_eq!(a.red_noise, b.red_noise);
        }
    }

    let event_grouped = Config {
        candidate_grouping: CandidateGrouping::Events,
        ..base
    };
    let event_serial = analyze(points.clone(), event_grouped).unwrap();
    for threads in [4, 24] {
        let event_parallel = analyze(
            points.clone(),
            Config {
                threads,
                ..event_grouped
            },
        )
        .unwrap();
        assert_eq!(
            render_csv(&event_serial.candidates),
            render_csv(&event_parallel.candidates)
        );
    }
}

#[test]
fn scoring_names_parse_and_invalid_options_are_rejected() {
    assert_eq!("white".parse::<Scoring>(), Ok(Scoring::White));
    assert_eq!("correlated".parse::<Scoring>(), Ok(Scoring::Correlated));
    assert!("pink".parse::<Scoring>().is_err());
    assert_eq!(
        "ephemeris".parse::<CandidateGrouping>(),
        Ok(CandidateGrouping::Ephemeris)
    );
    assert_eq!(
        "events".parse::<CandidateGrouping>(),
        Ok(CandidateGrouping::Events)
    );
    assert!("probability".parse::<CandidateGrouping>().is_err());
    let points = curve(3, false);
    for bad in [
        Config {
            min_transit_score: 0.0,
            ..Config::default()
        },
        Config {
            transit_candidates: 0,
            ..Config::default()
        },
    ] {
        assert!(analyze(points.clone(), bad).is_err());
    }
}

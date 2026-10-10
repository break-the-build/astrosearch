//! Light curves that combine TESS-like sectors with different cadences and gaps.

use astrosearch::{analyze, Candidate, Config, Point, Report};

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

const DURATION: f64 = 2.0 / 24.0;

/// A 30-minute sector at 250 ppm per sample, then a 200-second sector with equal
/// noise per unit time (750 ppm). Each has a 1-day orbit gap.
fn mixed(seed: u64, days_per_orbit: f64, errors: bool) -> Vec<Point> {
    let mut noise = Noise(seed);
    let mut points = Vec::new();
    for (start, seconds, sigma) in [(0.0, 1800.0, 250e-6), (40.0, 200.0, 750e-6)] {
        for orbit in [start, start + days_per_orbit + 1.0] {
            for i in 0..(days_per_orbit * 86400.0 / seconds) as usize {
                points.push(Point {
                    time: orbit + i as f64 * seconds / 86400.0,
                    flux: 1.0 + sigma * noise.gaussian(),
                    error: if errors { sigma } else { 0.0 },
                });
            }
        }
    }
    points
}

fn inject(points: &mut [Point], period: f64, center: f64, depth: f64) {
    for p in points {
        let phase = (p.time - center + period / 2.0).rem_euclid(period) - period / 2.0;
        if phase.abs() < DURATION / 2.0 {
            p.flux -= depth;
        }
    }
}

fn config(threads: usize) -> Config {
    Config {
        min_period: 2.0,
        max_period: 4.0,
        periods: 200,
        threads,
        ..Config::default()
    }
}

fn best_dip(report: &Report) -> Option<&Candidate> {
    report.candidates.iter().find(|c| c.kind.is_transit())
}

#[test]
fn sectors_report_their_own_cadence_and_noise() {
    for errors in [true, false] {
        let report = analyze(mixed(3, 6.0, errors), config(1)).unwrap();
        // Two sectors, each split by its orbit gap.
        assert_eq!(report.segments.len(), 4, "{:?}", report.segments);
        for (segment, (seconds, ppm)) in report.segments.iter().zip([
            (1800.0, 250.0),
            (1800.0, 250.0),
            (200.0, 750.0),
            (200.0, 750.0),
        ]) {
            assert!((segment.cadence_days * 86400.0 / seconds - 1.0).abs() < 1e-6);
            assert!(
                (segment.noise_ppm / ppm - 1.0).abs() < 0.2,
                "{errors}: {segment:?}"
            );
        }
    }
}

#[test]
fn shallow_transit_needs_the_low_cadence_sector_weighted_correctly() {
    // Both sectors carry equal transit signal-to-noise. Before per-segment noise,
    // the 30-minute sector was floored to 200-second scatter and lost ~9x weight.
    let (period, center, depth) = (2.73, 1.1, 210e-6);
    let mut points = mixed(11, 12.7, true);
    inject(&mut points, period, center, depth);
    let report = analyze(points, config(0)).unwrap();
    let best = best_dip(&report).expect("transit should be detected");
    assert!(
        (best.period - period).abs() * report.span_days / period < DURATION,
        "{best:?}"
    );
}

#[test]
fn mixed_cadence_noise_does_not_become_periodic_dips() {
    for seed in 0..4 {
        let report = analyze(mixed(500 + seed, 12.7, true), config(0)).unwrap();
        assert!(
            best_dip(&report).is_none(),
            "{seed}: {:?}",
            best_dip(&report)
        );
    }
}

#[test]
fn mixed_cadence_serial_and_parallel_scans_are_bit_identical() {
    let mut points = mixed(21, 6.0, false);
    inject(&mut points, 3.1, 0.4, 400e-6);
    for p in points
        .iter_mut()
        .filter(|p| (44.0..44.03).contains(&p.time))
    {
        p.flux += 0.01;
    }
    let serial = analyze(points.clone(), config(1)).unwrap();
    assert!(best_dip(&serial).is_some());
    for threads in [3, 24] {
        let parallel = analyze(points.clone(), config(threads)).unwrap();
        assert_eq!(parallel.warnings, serial.warnings);
        assert_eq!(
            astrosearch::render_csv(&parallel.candidates),
            astrosearch::render_csv(&serial.candidates)
        );
        let bits = |c: &Candidate| {
            [c.score, c.period, c.epoch, c.duration, c.depth, c.amplitude].map(f64::to_bits)
        };
        assert!(serial
            .candidates
            .iter()
            .zip(&parallel.candidates)
            .all(|(a, b)| bits(a) == bits(b)
                && a.odd_even_sigma.map(f64::to_bits) == b.odd_even_sigma.map(f64::to_bits)));
        for (a, b) in serial.segments.iter().zip(&parallel.segments) {
            assert_eq!(a.noise_ppm.to_bits(), b.noise_ppm.to_bits());
        }
    }
}

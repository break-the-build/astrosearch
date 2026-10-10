//! Mixed-cadence injection/recovery and null controls.
//!
//! Simulates two TESS-like sectors ~3 sectors apart: a 30-minute FFI sector
//! (250 ppm per sample) and a 200-second sector with equal noise per unit time
//! (750 ppm per sample), each with a mid-sector orbit gap. Box transits of known
//! period, phase, and depth are injected; null curves contain noise only. An
//! optional AR(1) component models time-correlated noise within each orbit.
//! Uses only the public API, so the same file can run against older builds.
//!
//!   cargo run --release --example mixed_cadence_injection -- [threads] [trials-per-depth] [null-trials] [red-rms-fraction]

use astrosearch::{analyze, Candidate, Config, Point};
use std::time::Instant;

const DURATION: f64 = 2.0 / 24.0;
const DEPTHS_PPM: [f64; 5] = [150.0, 200.0, 250.0, 300.0, 400.0];
const DEFAULT_NULL_TRIALS: u64 = 200;
const SCORE_THRESHOLDS: [f64; 8] = [6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 20.0, 25.0];

struct Rng(u64);

impl Rng {
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

/// (start day, cadence seconds, per-sample sigma) for each sector.
const SECTORS: [(f64, f64, f64); 2] = [(0.0, 1800.0, 250e-6), (90.0, 200.0, 750e-6)];

fn light_curve(seed: u64, red_fraction: f64) -> Vec<Point> {
    let mut rng = Rng(seed);
    let mut points = Vec::new();
    for (start, seconds, sigma) in SECTORS {
        // Two 13.2-day orbits separated by a 1-day downlink gap.
        for orbit_start in [start, start + 14.2] {
            let count = (13.2 * 86400.0 / seconds) as usize;
            let mut red = 0.0;
            for i in 0..count {
                red = 0.8 * red + 0.6 * rng.gaussian();
                let white = (1.0 - red_fraction * red_fraction).sqrt() * rng.gaussian();
                points.push(Point {
                    time: orbit_start + i as f64 * seconds / 86400.0,
                    flux: 1.0 + sigma * (white + red_fraction * red),
                    error: sigma,
                });
            }
        }
    }
    points
}

fn inject(points: &mut [Point], period: f64, center: f64, depth_ppm: f64) {
    for p in points {
        let phase = (p.time - center + period / 2.0).rem_euclid(period) - period / 2.0;
        if phase.abs() < DURATION / 2.0 {
            p.flux -= depth_ppm * 1e-6;
        }
    }
}

fn matches(c: &Candidate, period: f64, center: f64, span: f64) -> bool {
    let drift = (c.period - period).abs() * span / period;
    let phase = ((c.epoch - center + period / 2.0).rem_euclid(period) - period / 2.0).abs();
    drift <= DURATION && phase <= DURATION
}

fn dips(candidates: &[Candidate]) -> Vec<&Candidate> {
    candidates.iter().filter(|c| c.kind.is_transit()).collect()
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let argument = |i: usize, fallback: &str| {
        args.get(i)
            .map_or_else(|| fallback.to_string(), Clone::clone)
    };
    let threads: usize = argument(0, "0")
        .parse()
        .expect("threads must be an integer");
    let per_depth: u64 = argument(1, "12")
        .parse()
        .expect("trials-per-depth must be an integer");
    let null_default = DEFAULT_NULL_TRIALS.to_string();
    let null_trials: u64 = argument(2, &null_default)
        .parse()
        .expect("null-trials must be an integer");
    let red_fraction: f64 = argument(3, "0.5")
        .parse()
        .expect("red-rms-fraction must be numeric");
    assert!(
        per_depth > 0 && null_trials > 0 && (0.0..=1.0).contains(&red_fraction),
        "trial counts must be positive and red-rms-fraction in [0, 1]"
    );
    let config = Config {
        min_period: 1.0,
        max_period: 8.0,
        threads,
        ..Config::default()
    };
    let started = Instant::now();
    println!("# red_rms_fraction={red_fraction:.3}, ar1_rho=0.8");
    println!("depth_ppm,trials,top1_recovered,listed_recovered,median_matched_score");
    let mut matched_scores = Vec::new();
    for (d, depth) in DEPTHS_PPM.iter().enumerate() {
        let (mut top1, mut listed, mut scores) = (0, 0, Vec::new());
        for trial in 0..per_depth {
            let seed = 1000 * (d as u64 + 1) + trial;
            let mut rng = Rng(seed ^ 0x9E37_79B9_7F4A_7C15);
            let period = 1.5 + 4.5 * rng.uniform();
            let center = 0.3 + period * rng.uniform();
            let mut points = light_curve(seed, red_fraction);
            inject(&mut points, period, center, *depth);
            let report = analyze(points, config).expect("analysis");
            let found = dips(&report.candidates);
            if found
                .first()
                .is_some_and(|c| matches(c, period, center, report.span_days))
            {
                top1 += 1;
            }
            if let Some(c) = found
                .iter()
                .find(|c| matches(c, period, center, report.span_days))
            {
                listed += 1;
                scores.push(c.score);
                matched_scores.push((*depth, c.score));
            }
        }
        scores.sort_by(f64::total_cmp);
        let median = scores.get(scores.len() / 2).copied().unwrap_or(f64::NAN);
        println!("{depth},{per_depth},{top1},{listed},{median:.2}");
    }
    let mut null_max_scores = Vec::with_capacity(null_trials as usize);
    for trial in 0..null_trials {
        let report = analyze(light_curve(900_000 + trial, red_fraction), config).expect("analysis");
        let found = dips(&report.candidates);
        null_max_scores.push(found.iter().map(|c| c.score).fold(0.0_f64, f64::max));
    }
    println!("score_threshold,depth_ppm,recovered_injections,injections,null_alarm_curves,null_curves,recovery_fraction,false_alarm_fraction");
    for threshold in SCORE_THRESHOLDS {
        let false_alarm = null_max_scores
            .iter()
            .filter(|&&score| score >= threshold)
            .count();
        for depth in DEPTHS_PPM {
            let recovered = matched_scores
                .iter()
                .filter(|(d, score)| *d == depth && *score >= threshold)
                .count();
            let recovery = recovered as f64 / per_depth as f64;
            println!("{threshold:.1},{depth:.0},{recovered},{per_depth},{false_alarm},{null_trials},{recovery:.4},{:.4}",
                false_alarm as f64 / null_trials as f64);
        }
    }
    let any_alarm = null_max_scores
        .iter()
        .filter(|&&score| score >= 6.0)
        .count();
    let max_null_score = null_max_scores.iter().copied().fold(0.0_f64, f64::max);
    println!("# null_trials={null_trials}, score_ge_6_alarm_fraction={:.4}, max_null_score={max_null_score:.2}", any_alarm as f64 / null_trials as f64);
    println!("elapsed_seconds,{:.1}", started.elapsed().as_secs_f64());
}

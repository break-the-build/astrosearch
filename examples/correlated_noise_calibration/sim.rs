//! Simulated TESS-like light curves, real-curve loading, and transit injection.

use astrosearch::Point;
use std::{fs, path::Path};

pub struct Rng(pub u64);

impl Rng {
    pub fn uniform(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (((self.0 >> 11) as f64) + 0.5) / ((1_u64 << 53) as f64)
    }
    pub fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.uniform()
    }
    pub fn gaussian(&mut self) -> f64 {
        (-2.0 * self.uniform().ln()).sqrt() * (std::f64::consts::TAU * self.uniform()).cos()
    }
    /// Small Poisson draws by inversion.
    pub fn poisson(&mut self, mean: f64) -> usize {
        let (limit, mut product, mut count) = ((-mean).exp(), self.uniform(), 0);
        while product > limit {
            product *= self.uniform();
            count += 1;
        }
        count
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum NoiseModel {
    White,
    Correlated,
}

impl NoiseModel {
    pub fn name(self) -> &'static str {
        match self {
            Self::White => "white",
            Self::Correlated => "correlated",
        }
    }
}

/// A sector observed with one cadence: start day, cadence seconds, white sigma,
/// and the correlated-model OU amplitude/timescale and glitch amplitude range.
struct Sector {
    start: f64,
    seconds: f64,
    white: f64,
    red: f64,
    tau: f64,
    glitch: (f64, f64),
}

/// A 30-minute FFI sector, then a 10-minute FFI sector 90 days later, with
/// equal white noise per unit time. In the correlated model the 10-minute
/// sector carries most of the red noise, so a single global beta would be wrong.
const SECTORS: [Sector; 2] = [
    Sector {
        start: 0.0,
        seconds: 1800.0,
        white: 250e-6,
        red: 60e-6,
        tau: 0.08,
        glitch: (60e-6, 150e-6),
    },
    Sector {
        start: 90.0,
        seconds: 600.0,
        white: 433e-6,
        red: 160e-6,
        tau: 0.12,
        glitch: (120e-6, 300e-6),
    },
];
const ORBIT_DAYS: f64 = 13.2;
const ORBIT_GAP: f64 = 1.0;

/// Two TESS-like sectors, each two orbits. The correlated model adds an
/// Ornstein-Uhlenbeck process, scattered-light ramps at orbit starts, and
/// Gaussian dip-like glitches (about 1.5 per orbit), none of which the
/// quoted errors include.
pub fn simulate(model: NoiseModel, seed: u64) -> Vec<Point> {
    let mut rng = Rng(seed);
    let mut points = Vec::new();
    for sector in &SECTORS {
        for orbit in 0..2 {
            let start = sector.start + orbit as f64 * (ORBIT_DAYS + ORBIT_GAP);
            let count = (ORBIT_DAYS * 86400.0 / sector.seconds) as usize;
            let step = sector.seconds / 86400.0;
            let glitches: Vec<(f64, f64, f64)> = match model {
                NoiseModel::White => Vec::new(),
                NoiseModel::Correlated => (0..rng.poisson(1.5))
                    .map(|_| {
                        (
                            rng.range(start, start + ORBIT_DAYS),
                            rng.range(0.02, 0.06),
                            rng.range(sector.glitch.0, sector.glitch.1),
                        )
                    })
                    .collect(),
            };
            let ramp = rng.range(200e-6, 800e-6);
            let decay = (-step / sector.tau).exp();
            let mut state = 0.0;
            for i in 0..count {
                let time = start + i as f64 * step;
                let mut flux = 1.0 + sector.white * rng.gaussian();
                if model == NoiseModel::Correlated {
                    state =
                        state * decay + sector.red * (1.0 - decay * decay).sqrt() * rng.gaussian();
                    flux += state + ramp * (-(time - start) / 0.3).exp();
                    for &(center, width, depth) in &glitches {
                        flux -= depth * (-0.5 * ((time - center) / width).powi(2)).exp();
                    }
                }
                points.push(Point {
                    time,
                    flux,
                    error: sector.white,
                });
            }
        }
    }
    points
}

/// Reads `time,flux,flux_err` CSVs written by scripts/export_tess.py.
pub fn read_curve(path: &Path) -> Vec<Point> {
    let text = fs::read_to_string(path).expect("readable light curve");
    text.lines()
        .skip(1)
        .filter_map(|line| {
            let mut fields = line.split(',').map(|v| v.trim().parse::<f64>());
            match (fields.next(), fields.next(), fields.next()) {
                (Some(Ok(time)), Some(Ok(flux)), Some(Ok(error))) => {
                    Some(Point { time, flux, error })
                }
                _ => None,
            }
        })
        .collect()
}

/// Reflect normalized flux about its median: dips become brightenings.
pub fn invert(points: &[Point]) -> Vec<Point> {
    let mut fluxes: Vec<f64> = points.iter().map(|p| p.flux).collect();
    fluxes.sort_by(f64::total_cmp);
    let center = fluxes[fluxes.len() / 2];
    points
        .iter()
        .map(|p| Point {
            flux: 2.0 * center - p.flux,
            ..*p
        })
        .collect()
}

#[derive(Clone, Copy, Debug)]
pub struct Planet {
    pub period: f64,
    pub epoch: f64,
    pub duration: f64,
    pub depth: f64,
}

impl Planet {
    /// Signed time from the nearest transit center.
    pub fn offset(&self, time: f64) -> f64 {
        (time - self.epoch + self.period / 2.0).rem_euclid(self.period) - self.period / 2.0
    }

    /// Trapezoid with ingress and egress each 10% of the duration.
    fn shape(&self, time: f64) -> f64 {
        let x = self.offset(time).abs();
        let (half, ramp) = (self.duration / 2.0, 0.1 * self.duration);
        if x >= half {
            0.0
        } else if x <= half - ramp {
            1.0
        } else {
            (half - x) / ramp
        }
    }

    pub fn inject(&self, points: &mut [Point]) {
        for p in points {
            p.flux *= 1.0 - self.depth * self.shape(p.time);
        }
    }

    /// Transits with at least one sample inside them.
    pub fn observed_events(&self, points: &[Point]) -> usize {
        let mut cycles: Vec<i64> = points
            .iter()
            .filter(|p| self.offset(p.time).abs() < self.duration / 2.0)
            .map(|p| ((p.time - self.epoch) / self.period).round() as i64)
            .collect();
        cycles.dedup();
        cycles.len()
    }
}

/// Per-sample white sigma from neighbor differences within 3 cadences.
pub fn sample_sigma(points: &[Point]) -> f64 {
    let mut spacing: Vec<f64> = points.windows(2).map(|p| p[1].time - p[0].time).collect();
    spacing.sort_by(f64::total_cmp);
    let cadence = spacing[spacing.len() / 2];
    let mut differences: Vec<f64> = points
        .windows(2)
        .filter(|p| p[1].time - p[0].time <= 3.0 * cadence)
        .map(|p| p[1].flux - p[0].flux)
        .collect();
    differences.sort_by(f64::total_cmp);
    let median = differences[differences.len() / 2];
    let mut deviations: Vec<f64> = differences.iter().map(|d| (d - median).abs()).collect();
    deviations.sort_by(f64::total_cmp);
    1.4826 * deviations[deviations.len() / 2] / std::f64::consts::SQRT_2
}

/// Depth giving a target white-noise box SNR: SNR = depth * sqrt(sum 1/sigma^2)
/// over in-transit samples, with sigma the larger of quoted error and `floor`.
pub fn depth_for_snr(points: &[Point], planet: &Planet, snr: f64, floor: f64) -> f64 {
    let weight: f64 = points
        .iter()
        .filter(|p| planet.offset(p.time).abs() < planet.duration / 2.0)
        .map(|p| 1.0 / p.error.max(floor).powi(2))
        .sum();
    if weight > 0.0 {
        snr / weight.sqrt()
    } else {
        f64::NAN
    }
}

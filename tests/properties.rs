//! Randomized property tests for invariants the pipeline guarantees by
//! construction (no external crates; a fixed-seed generator keeps them
//! reproducible). Each property runs over several random light curves.

use astrosearch::{analyze, render_csv, render_csv_with_tier, Config, Point};

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
    fn below(&mut self, n: usize) -> usize {
        (self.uniform() * n as f64) as usize % n
    }
}

/// A light curve with gaps, mixed errors and an injected transit.
fn random_curve(rng: &mut Rng) -> Vec<Point> {
    let count = 900 + rng.below(600);
    let period = 2.0 + rng.uniform() * 2.0;
    let center = rng.uniform() * period;
    let depth = 0.001 + rng.uniform() * 0.003;
    let mut time = 0.0;
    (0..count)
        .map(|_| {
            time += if rng.uniform() < 0.01 { 1.5 } else { 0.02 };
            let phase = ((time - center + period / 2.0).rem_euclid(period) - period / 2.0).abs();
            let dip = if phase < 0.05 { depth } else { 0.0 };
            Point {
                time,
                flux: 1.0 + 0.0003 * rng.gaussian() - dip,
                error: 0.0002 + 0.0002 * rng.uniform(),
            }
        })
        .collect()
}

fn config() -> Config {
    Config {
        min_period: 1.5,
        max_period: 4.5,
        periods: 150,
        durations: 4,
        threads: 1,
        transit_search_passes: 2,
        ..Config::default()
    }
}

#[test]
fn candidates_do_not_depend_on_input_row_order() {
    let mut rng = Rng(11);
    for _ in 0..4 {
        let points = random_curve(&mut rng);
        let mut shuffled = points.clone();
        for i in (1..shuffled.len()).rev() {
            shuffled.swap(i, rng.below(i + 1));
        }
        let a = analyze(points, config()).unwrap();
        let b = analyze(shuffled, config()).unwrap();
        assert_eq!(render_csv(&a.candidates), render_csv(&b.candidates));
    }
}

#[test]
fn invalid_rows_never_change_candidates() {
    let mut rng = Rng(23);
    for _ in 0..4 {
        let points = random_curve(&mut rng);
        let mut polluted = points.clone();
        let span = points.last().unwrap().time;
        for _ in 0..50 {
            let time = rng.uniform() * span + 0.0013; // off the sampling grid
            let flux = [f64::NAN, f64::INFINITY, -1.0, 0.0][rng.below(4)];
            polluted.push(Point {
                time,
                flux,
                error: 0.0002,
            });
            polluted.push(Point {
                time: f64::NAN,
                flux: 1.0,
                error: 0.0002,
            });
        }
        let a = analyze(points, config()).unwrap();
        let b = analyze(polluted, config()).unwrap();
        assert_eq!(render_csv(&a.candidates), render_csv(&b.candidates));
    }
}

#[test]
fn csv_rows_have_header_width_and_finite_numbers() {
    let mut rng = Rng(37);
    for _ in 0..3 {
        let report = analyze(random_curve(&mut rng), config()).unwrap();
        let csv = render_csv_with_tier(&report.candidates, 10.0);
        let mut lines = csv.lines();
        let width = lines.next().unwrap().split(',').count();
        for line in lines {
            // Quoted cells contain no commas in numeric columns; split outside quotes.
            let mut cells = Vec::new();
            let (mut cell, mut quoted) = (String::new(), false);
            for ch in line.chars() {
                match ch {
                    '"' => quoted = !quoted,
                    ',' if !quoted => cells.push(std::mem::take(&mut cell)),
                    _ => cell.push(ch),
                }
            }
            cells.push(cell);
            assert_eq!(cells.len(), width, "{line}");
            for index in [3, 4, 5, 6, 7, 11] {
                let value: f64 = cells[index].parse().unwrap();
                assert!(value.is_finite(), "column {index} in {line}");
            }
        }
    }
}

use astrosearch::{analyze, render_csv, render_csv_with_tier, Config, Point, Report};

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

fn curve(count: usize, seed: u64) -> Vec<Point> {
    let mut noise = Noise(seed);
    (0..count)
        .map(|i| Point {
            time: i as f64 * 0.02,
            flux: 1.0 + 0.0002 * noise.gaussian(),
            error: 0.0002,
        })
        .collect()
}

fn inject(points: &mut [Point], period: f64, center: f64, duration: f64, depth: f64) {
    for point in points {
        let distance =
            ((point.time - center + period / 2.0).rem_euclid(period) - period / 2.0).abs();
        if distance < duration / 2.0 {
            point.flux -= depth;
        }
    }
}

fn config(min_period: f64, max_period: f64) -> Config {
    Config {
        min_period,
        max_period,
        periods: 200,
        durations: 6,
        threads: 1,
        ..Config::default()
    }
}

fn has_period(report: &Report, period: f64, tolerance: f64) -> bool {
    report
        .candidates
        .iter()
        .any(|c| c.kind.is_transit() && (c.period - period).abs() < tolerance)
}

#[test]
fn invalid_parameters_and_duplicate_only_data_return_errors() {
    for bad in [f64::NAN, f64::INFINITY, 0.0, -1.0, 1e-320] {
        assert!(analyze(
            curve(40, 1),
            Config {
                trend_days: bad,
                ..Config::default()
            }
        )
        .is_err());
    }
    assert!(analyze(
        curve(40, 1),
        Config {
            periods: usize::MAX,
            ..Config::default()
        }
    )
    .is_err());
    let repeated = vec![
        Point {
            time: 1.0,
            flux: 1.0,
            error: 0.1
        };
        100
    ];
    assert!(analyze(repeated, Config::default())
        .unwrap_err()
        .contains("distinct"));
}

#[test]
#[cfg(unix)]
fn cli_refuses_to_overwrite_input_or_its_hard_link() {
    let input = std::env::temp_dir().join(format!("astrosearch-review-{}.csv", std::process::id()));
    let linked = input.with_extension("linked.csv");
    let contents = "time,flux\n0,1\n1,1\n";
    std::fs::write(&input, contents).unwrap();
    for output in [&input, &linked] {
        if output == &linked {
            std::fs::hard_link(&input, &linked).unwrap();
        }
        let result = std::process::Command::new(env!("CARGO_BIN_EXE_astrosearch"))
            .arg("--input")
            .arg(&input)
            .arg("--output")
            .arg(output)
            .output()
            .unwrap();
        assert_eq!(result.status.code(), Some(2));
        assert!(String::from_utf8_lossy(&result.stderr).contains("refusing to overwrite"));
        assert_eq!(std::fs::read_to_string(&input).unwrap(), contents);
    }
    std::fs::remove_file(&linked).unwrap();
    std::fs::remove_file(&input).unwrap();
}

#[test]
fn cli_refuses_candidate_and_periodogram_path_collision() {
    let directory =
        std::env::temp_dir().join(format!("astrosearch-collision-{}", std::process::id()));
    std::fs::create_dir_all(&directory).unwrap();
    let input = directory.join("input.csv");
    let output = directory.join("result.csv");
    let alias = directory.join(".").join("result.csv");
    std::fs::write(&input, "time,flux\n0,1\n1,1\n").unwrap();
    let result = std::process::Command::new(env!("CARGO_BIN_EXE_astrosearch"))
        .arg("--input")
        .arg(&input)
        .arg("--output")
        .arg(&output)
        .arg("--periodogram")
        .arg(&alias)
        .output()
        .unwrap();
    assert_eq!(result.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&result.stderr)
        .contains("candidate output and periodogram output must be different files"));
    std::fs::remove_dir_all(directory).unwrap();
}

#[test]
fn weak_repeated_transits_accumulate_significance() {
    let mut points = curve(1500, 17);
    inject(&mut points, 3.0, 0.37, 0.12, 0.00023);
    let mut sums = [0.0; 10];
    let mut counts = [0_usize; 10];
    for point in &points {
        if (point.time.rem_euclid(3.0) - 0.37).abs() < 0.06 {
            let cycle = (point.time / 3.0).floor() as usize;
            sums[cycle] += 1.0 - point.flux;
            counts[cycle] += 1;
        }
    }
    let individually_significant = sums
        .iter()
        .zip(counts)
        .filter(|(sum, count)| *count > 0 && **sum / (0.0002 * (*count as f64).sqrt()) >= 3.0)
        .count();
    // This fixture would fail the previous "80% of epochs >= 3 sigma" gate.
    assert!(individually_significant * 5 < counts.len() * 4);
    let report = analyze(points, config(2.7, 3.3)).unwrap();
    assert!(has_period(&report, 3.0, 0.02), "{:?}", report.candidates);
}

#[test]
fn nearby_planet_periods_are_both_recovered() {
    let mut points = curve(3000, 83);
    inject(&mut points, 3.0, 0.37, 0.12, 0.004);
    inject(&mut points, 3.3, 1.13, 0.12, 0.003);
    let report = analyze(points, config(2.8, 3.6)).unwrap();
    assert!(has_period(&report, 3.0, 0.015), "{:?}", report.candidates);
    assert!(has_period(&report, 3.3, 0.015), "{:?}", report.candidates);
}

#[test]
fn two_transits_can_have_period_longer_than_half_baseline() {
    let mut points = curve(1350, 31);
    inject(&mut points, 18.0, 2.07, 0.14, 0.006);
    let report = analyze(points, config(17.0, 19.0)).unwrap();
    assert!(has_period(&report, 18.0, 0.04), "{:?}", report.candidates);
    let best = report
        .candidates
        .iter()
        .find(|c| c.kind.is_transit())
        .unwrap();
    assert_eq!(best.n_events, 2);
    // Two observed transits are reported as a duo, not a periodic candidate.
    assert_eq!(best.kind, astrosearch::CandidateKind::DuoTransit);
    assert!(best.note.contains("period ambiguous"));
    assert!((best.epoch - 2.07).abs() < 0.08);
}

#[test]
fn trend_method_none_only_removes_segment_levels() {
    // Pre-flattened data with a sector offset: no detrending, still found.
    let mut points = curve(1500, 47);
    for p in points.iter_mut().filter(|p| p.time > 15.0) {
        p.time += 5.0;
        p.flux *= 1.02;
    }
    inject(&mut points, 3.1, 0.6, 0.12, 0.003);
    let cfg = Config {
        trend_method: astrosearch::TrendMethod::None,
        ..config(2.8, 3.4)
    };
    assert!(has_period(&analyze(points, cfg).unwrap(), 3.1, 0.03));
    assert_eq!("none".parse(), Ok(astrosearch::TrendMethod::None));
}

#[test]
fn flux_scale_and_absolute_time_origin_do_not_change_discovery() {
    let mut points = curve(1200, 27);
    inject(&mut points, 3.1, 0.37, 0.12, 0.004);
    let original = analyze(points.clone(), config(2.9, 3.3)).unwrap();
    for point in &mut points {
        point.flux *= 1000.0;
        point.error *= 1000.0;
        point.time += 2_457_000.0;
    }
    let shifted = analyze(points, config(2.9, 3.3)).unwrap();
    let a = original
        .candidates
        .iter()
        .find(|c| c.kind.is_transit())
        .unwrap();
    let b = shifted
        .candidates
        .iter()
        .find(|c| c.kind.is_transit())
        .unwrap();
    assert!((a.period - b.period).abs() < 0.003);
    assert!((a.score / b.score - 1.0).abs() < 0.05);
    assert!((a.epoch - (b.epoch - 2_457_000.0)).abs() < 0.06);
}

#[test]
fn null_curve_and_single_events_do_not_become_periodic_planets() {
    let points = curve(900, 8080);
    let cfg = config(1.5, 4.5);
    let null = analyze(points.clone(), cfg).unwrap();
    assert!(null.candidates.is_empty(), "{:?}", null.candidates);
    for depth in [-0.04, 0.008] {
        let mut event_points = points.clone();
        for point in &mut event_points {
            if (10.0..10.06).contains(&point.time) {
                point.flux -= depth;
            }
        }
        let report = analyze(event_points, cfg).unwrap();
        assert!(
            !report.candidates.iter().any(|c| c.kind.is_transit()),
            "{:?}",
            report.candidates
        );
        let kind = if depth > 0.0 {
            "isolated_dip"
        } else {
            "flare_or_brightening"
        };
        assert!(
            report.candidates.iter().any(|c| c.kind == kind),
            "{:?}",
            report.candidates
        );
    }
}

#[test]
fn variability_respects_bounds_and_recovers_an_allowed_sinusoid() {
    let mut points = curve(1200, 90);
    for point in &mut points {
        point.flux += 0.004 * (std::f64::consts::TAU * point.time / 2.4).sin();
    }
    let bounded = analyze(points.clone(), config(2.6, 3.4)).unwrap();
    assert!(bounded
        .candidates
        .iter()
        .filter(|c| c.kind == "periodic_variability")
        .all(|c| c.period >= 2.6 - 1e-12 && c.period <= 3.4 + 1e-12));
    let allowed = analyze(points, config(1.9, 2.8)).unwrap();
    assert!(
        allowed
            .candidates
            .iter()
            .any(|c| c.kind == "periodic_variability" && (c.period - 2.4).abs() < 0.02),
        "{:?}",
        allowed.candidates
    );
}

#[test]
fn serial_and_parallel_scans_return_identical_candidates() {
    let mut points = curve(1200, 44);
    inject(&mut points, 3.0, 0.37, 0.12, 0.008);
    for point in &mut points {
        point.flux += 0.003 * (std::f64::consts::TAU * point.time / 3.3).sin();
        if (10.0..10.08).contains(&point.time) {
            point.flux += 0.02;
        }
    }
    points.retain(|p| !(7.0..8.0).contains(&p.time));
    let cfg = config(2.8, 3.6);
    let serial = analyze(points.clone(), cfg).unwrap();
    assert!(serial.candidates.iter().any(|c| c.kind.is_transit()));
    assert!(serial
        .candidates
        .iter()
        .any(|c| c.kind == "periodic_variability"));
    for threads in [4, 24] {
        let parallel = analyze(points.clone(), Config { threads, ..cfg }).unwrap();
        assert_eq!(parallel.warnings, serial.warnings);
        assert_eq!(
            astrosearch::render_csv(&parallel.candidates),
            astrosearch::render_csv(&serial.candidates)
        );
        for (a, b) in serial.candidates.iter().zip(&parallel.candidates) {
            let numbers = |c: &astrosearch::Candidate| {
                [c.score, c.period, c.epoch, c.duration, c.depth, c.amplitude].map(f64::to_bits)
            };
            assert_eq!(numbers(a), numbers(b));
            assert_eq!(
                a.odd_even_sigma.map(f64::to_bits),
                b.odd_even_sigma.map(f64::to_bits)
            );
        }
    }
}

#[test]
fn masked_transits_are_excluded_not_counted_as_flat_data() {
    // Planet B (3 d) shares every other transit time with planet A (2 d).
    let base = curve(3000, 83);
    let mut both = base.clone();
    inject(&mut both, 2.0, 0.5, 0.12, 0.004);
    inject(&mut both, 3.0, 0.5, 0.12, 0.0012);
    let mut alone = base;
    inject(&mut alone, 3.0, 0.5, 0.12, 0.0012);
    let cfg = Config {
        scoring: astrosearch::Scoring::White,
        transit_search_passes: 2,
        ..config(1.5, 4.0)
    };
    let score_near = |report: &Report, period: f64| {
        report
            .candidates
            .iter()
            .filter(|c| c.kind.is_transit() && (c.period - period).abs() < 0.01)
            .map(|c| c.white_score)
            .fold(f64::NAN, f64::max)
    };
    let solo = score_near(&analyze(alone, cfg).unwrap(), 3.0);
    let masked = score_near(&analyze(both, cfg).unwrap(), 3.0);
    // Excluding A's windows leaves B half its events: about 1/sqrt(2) of its
    // solo score. Replacing them with flat data would halve it instead.
    let ratio = masked / solo;
    assert!(
        (0.6..0.85).contains(&ratio),
        "masked {masked:.1} solo {solo:.1}"
    );
}

#[test]
fn candidate_limit_not_number_of_search_passes_bounds_transit_output() {
    let points = curve(1500, 71);
    let dips = |passes: usize, limit: usize| {
        let cfg = Config {
            min_transit_score: 3.0,
            transit_search_passes: passes,
            transit_candidates: limit,
            ..config(1.0, 8.0)
        };
        analyze(points.clone(), cfg)
            .unwrap()
            .candidates
            .iter()
            .filter(|c| c.kind.is_transit())
            .count()
    };
    // Low-threshold noise windows are plentiful, so the slots should fill.
    assert!(dips(2, 12) > 2, "output stopped at the number of passes");
    assert!(dips(2, 12) <= 12);
    assert!(dips(5, 4) <= 4);
}

#[test]
fn tier_column_separates_detection_from_screening_transits() {
    let mut points = curve(3000, 83);
    inject(&mut points, 3.0, 0.37, 0.12, 0.004);
    let report = analyze(points, config(2.8, 3.6)).unwrap();
    let best = report
        .candidates
        .iter()
        .filter(|c| c.kind.is_transit())
        .map(|c| c.jackknife_score())
        .fold(0.0, f64::max);
    assert!(best > 0.0, "{:?}", report.candidates);
    let csv = render_csv_with_tier(&report.candidates, best);
    let mut lines = csv.lines();
    assert!(lines.next().unwrap().ends_with(",jackknife_score,tier"));
    let tiers: Vec<&str> = lines.map(|l| l.rsplit(',').next().unwrap()).collect();
    assert_eq!(tiers.iter().filter(|t| **t == "detection").count(), 1);
    // A threshold above every score leaves only screening-tier transits.
    let none = render_csv_with_tier(&report.candidates, best + 1.0);
    assert!(!none.contains(",detection"));
    // The plain renderer keeps the original header for existing consumers.
    assert!(render_csv(&report.candidates)
        .lines()
        .next()
        .unwrap()
        .ends_with(",note"));
}

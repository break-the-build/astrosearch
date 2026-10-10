//! Matched false-alarm comparison of white versus correlated transit scoring.
//!
//! Null sets (white and correlated simulated noise, inverted real TESS curves)
//! set each statistic's threshold at a fixed false-alarm rate, then injection
//! families are scored against those thresholds: varied durations, two-event
//! planets, two-planet systems, and injections into real curves.
//!
//!   cargo run --release --example correlated_noise_calibration -- \
//!       [--scale F] [--real DIR] [--out DIR]

mod sim;

use astrosearch::{
    analyze, render_csv, Candidate, CandidateGrouping, Config, Point, Report, Scoring,
};
use sim::{depth_for_snr, invert, read_curve, sample_sigma, simulate, NoiseModel, Planet, Rng};
use std::{fmt::Write as _, fs, path::PathBuf, time::Instant};

const MODES: [Scoring; 2] = [Scoring::White, Scoring::Correlated];
const ALPHAS: [f64; 3] = [0.10, 0.05, 0.01];
const MODELS: [NoiseModel; 2] = [NoiseModel::White, NoiseModel::Correlated];

fn config(mode: Scoring, min_period: f64, max_period: f64, grouping: CandidateGrouping) -> Config {
    Config {
        min_period,
        max_period,
        scoring: mode,
        candidate_grouping: grouping,
        // Low enough that null curves still report their strongest dip.
        min_transit_score: 4.0,
        ..Config::default()
    }
}

/// One curve searched by both statistics.
struct Outcome {
    top: [f64; 2],
    reports: [Report; 2],
}

struct Harness {
    seconds: [f64; 2],
    raw: String,
    determinism_checks: usize,
    determinism_failures: usize,
    grouping: CandidateGrouping,
}

impl Harness {
    fn search(&mut self, points: &[Point], range: (f64, f64), check: bool) -> Outcome {
        let reports = MODES.map(|mode| {
            let started = Instant::now();
            let report = analyze(
                points.to_vec(),
                config(mode, range.0, range.1, self.grouping),
            )
            .expect("analysis");
            self.seconds[mode as usize] += started.elapsed().as_secs_f64();
            report
        });
        if check {
            let serial = Config {
                threads: 1,
                ..config(Scoring::Correlated, range.0, range.1, self.grouping)
            };
            let report = analyze(points.to_vec(), serial).expect("analysis");
            self.determinism_checks += 1;
            if render_csv(&report.candidates) != render_csv(&reports[1].candidates) {
                self.determinism_failures += 1;
            }
        }
        let top = [0, 1].map(|m| dips(&reports[m]).first().map_or(0.0, |c| c.score));
        Outcome { top, reports }
    }
}

fn dips(report: &Report) -> Vec<&Candidate> {
    report
        .candidates
        .iter()
        .filter(|c| c.kind.is_transit())
        .collect()
}

/// A candidate ephemeris that predicts every listed transit time within a duration.
fn predicts_ephemeris(period: f64, epoch: f64, times: &[f64], duration: f64) -> bool {
    times
        .iter()
        .all(|t| ((t - epoch + period / 2.0).rem_euclid(period) - period / 2.0).abs() <= duration)
}

fn ephemeris_matches_planet(period: f64, epoch: f64, planet: &Planet, span: f64) -> bool {
    (period - planet.period).abs() * span / planet.period <= planet.duration
        && predicts_ephemeris(period, epoch, &[planet.epoch], planet.duration)
}

fn ephemerides(c: &Candidate) -> impl Iterator<Item = (f64, f64, f64, f64)> + '_ {
    std::iter::once((c.period, c.epoch, c.duration, c.score)).chain(
        c.alias_ephemerides
            .iter()
            .map(|a| (a.period, a.epoch, a.duration, a.score)),
    )
}

fn matched_ephemeris_score(c: &Candidate, planet: &Planet, span: f64) -> Option<f64> {
    ephemerides(c)
        .find(|(period, epoch, _, _)| ephemeris_matches_planet(*period, *epoch, planet, span))
        .map(|(_, _, _, score)| score)
}

/// Thresholds where at most `alpha` of nulls score strictly above.
fn thresholds(nulls: &[f64]) -> Vec<(f64, f64, f64)> {
    let mut sorted = nulls.to_vec();
    sorted.sort_by(|a, b| b.total_cmp(a));
    ALPHAS
        .iter()
        .map(|&alpha| {
            let k = (alpha * sorted.len() as f64).floor() as usize;
            let t = sorted[k.min(sorted.len() - 1)];
            let achieved = sorted.iter().filter(|&&s| s > t).count() as f64 / sorted.len() as f64;
            (alpha, t, achieved)
        })
        .collect()
}

fn wilson(successes: usize, trials: usize) -> (f64, f64) {
    if trials == 0 {
        return (0.0, 1.0);
    }
    let z = 1.96;
    let n = trials as f64;
    let p = successes as f64 / n;
    let denominator = 1.0 + z * z / n;
    let center = (p + z * z / (2.0 * n)) / denominator;
    let half = z * ((p * (1.0 - p) / n + z * z / (4.0 * n * n)).sqrt()) / denominator;
    ((center - half).max(0.0), (center + half).min(1.0))
}

/// A detection: the matched candidate's score must exceed the null threshold.
#[derive(Clone, Copy)]
struct Trial {
    /// Rank of the matching candidate group and its ephemeris score, per mode.
    matched: [Option<Matched>; 2],
}

#[derive(Clone, Copy)]
struct Matched {
    score: f64,
    rank: usize,
}

fn recovered(trials: &[Trial], mode: usize, threshold: f64, max_rank: usize) -> usize {
    trials
        .iter()
        .filter(|t| t.matched[mode].is_some_and(|m| m.score > threshold && m.rank <= max_rank))
        .count()
}

fn table_row(out: &mut String, label: &str, trials: &[Trial], cuts: &[Vec<(f64, f64, f64)>; 2]) {
    let mut cells = Vec::new();
    for (a, _) in ALPHAS.iter().enumerate() {
        let w = recovered(trials, 0, cuts[0][a].1, 12);
        let r = recovered(trials, 1, cuts[1][a].1, 12);
        let (wl, wh) = wilson(w, trials.len());
        let (rl, rh) = wilson(r, trials.len());
        cells.push(format!(
            "{w}/{} [{wl:.2}, {wh:.2}] → **{r}/{} [{rl:.2}, {rh:.2}]**",
            trials.len(),
            trials.len()
        ));
    }
    let ranks = [1, 5, 12].map(|rank| {
        let w = recovered(trials, 0, cuts[0][1].1, rank);
        let r = recovered(trials, 1, cuts[1][1].1, rank);
        let (wl, wh) = wilson(w, trials.len());
        let (rl, rh) = wilson(r, trials.len());
        format!(
            "{w}/{} [{wl:.2}, {wh:.2}] → **{r}/{} [{rl:.2}, {rh:.2}]**",
            trials.len(),
            trials.len()
        )
    });
    writeln!(
        out,
        "| {label} | {} | {} | {} | {} | {} |",
        trials.len(),
        cells.join(" | "),
        ranks[0],
        ranks[1],
        ranks[2]
    )
    .unwrap();
}

fn threshold_rows(out: &mut String, name: &str, nulls: &[[f64; 2]]) -> [Vec<(f64, f64, f64)>; 2] {
    let cuts = [0, 1].map(|m| thresholds(&nulls.iter().map(|n| n[m]).collect::<Vec<_>>()));
    let fixed = [0, 1].map(|m| nulls.iter().filter(|n| n[m] >= 6.0).count());
    let cells: Vec<String> = (0..ALPHAS.len())
        .map(|a| format!("{:.2} / {:.2}", cuts[0][a].1, cuts[1][a].1))
        .collect();
    let interval = |k: usize| {
        let (lo, hi) = wilson(k, nulls.len());
        format!("{k} [{lo:.3}, {hi:.3}]")
    };
    writeln!(
        out,
        "| {name} | {} | {} | {} / {} |",
        nulls.len(),
        cells.join(" | "),
        interval(fixed[0]),
        interval(fixed[1])
    )
    .unwrap();
    cuts
}

fn single(outcome: &Outcome, accept: impl Fn(f64, f64, f64) -> bool) -> Trial {
    Trial {
        matched: [0, 1].map(|m| {
            dips(&outcome.reports[m])
                .iter()
                .enumerate()
                .find_map(|(index, candidate)| {
                    ephemerides(candidate)
                        .find(|(period, epoch, duration, _)| accept(*period, *epoch, *duration))
                        .map(|(_, _, _, score)| Matched {
                            score,
                            rank: index + 1,
                        })
                })
        }),
    }
}

fn record(
    h: &mut Harness,
    family: &str,
    model: &str,
    id: usize,
    detail: &str,
    o: &Outcome,
    t: &[Trial],
) {
    let fmt =
        |v: Option<Matched>| v.map_or(String::new(), |m| format!("{:.4}@{}", m.score, m.rank));
    let matched: Vec<String> = t.iter().flat_map(|t| t.matched.map(fmt)).collect();
    writeln!(
        h.raw,
        "{family},{model},{id},{detail},{:.4},{:.4},{}",
        o.top[0],
        o.top[1],
        matched.join(";")
    )
    .unwrap();
}

fn main() {
    let mut scale = 1.0;
    let mut real_dir = PathBuf::from("data/calibration/s32_tess_spoc_600");
    let mut real_limit = usize::MAX;
    let mut out_dir = PathBuf::from("target/calibration");
    let mut grouping = CandidateGrouping::Ephemeris;
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        let value = args.next().expect("option value");
        match arg.as_str() {
            "--scale" => scale = value.parse().expect("numeric scale"),
            "--real" => real_dir = PathBuf::from(value),
            "--real-limit" => real_limit = value.parse().expect("integer real-curve limit"),
            "--out" => out_dir = PathBuf::from(value),
            "--grouping" => grouping = value.parse().expect("ephemeris or events grouping"),
            other => panic!("unknown option {other}"),
        }
    }
    let n = |count: f64| ((count * scale).round() as usize).max(4);
    let started = Instant::now();
    let mut h = Harness {
        seconds: [0.0; 2],
        raw: String::from(
            "family,noise,id,detail,top_white,top_correlated,matched_white;matched_correlated...\n",
        ),
        determinism_checks: 0,
        determinism_failures: 0,
        grouping,
    };
    let mut report = String::new();
    let mut nulls_short = Vec::new();
    let mut nulls_long = Vec::new();
    writeln!(
        report,
        "### Candidate grouping: {:?}\n\n### Null thresholds (white / correlated statistic)\n",
        grouping
    )
    .unwrap();
    writeln!(
        report,
        "| Null set | Curves | T at 10% | T at 5% | T at 1% | Curves ≥ 6.0 |"
    )
    .unwrap();
    writeln!(report, "| --- | ---: | ---: | ---: | ---: | ---: |").unwrap();
    for (m, model) in MODELS.iter().enumerate() {
        let short: Vec<[f64; 2]> = (0..n(300.0))
            .map(|i| {
                let o = h.search(
                    &simulate(*model, 10_000 * (m as u64 + 1) + i as u64),
                    (1.0, 10.0),
                    i < 3,
                );
                record(&mut h, "null_1-10d", model.name(), i, "", &o, &[]);
                o.top
            })
            .collect();
        let long: Vec<[f64; 2]> = (0..n(150.0))
            .map(|i| {
                let o = h.search(
                    &simulate(*model, 50_000 * (m as u64 + 1) + i as u64),
                    (20.0, 50.0),
                    false,
                );
                record(&mut h, "null_20-50d", model.name(), i, "", &o, &[]);
                o.top
            })
            .collect();
        nulls_short.push(threshold_rows(
            &mut report,
            &format!("{} noise, P 1–10 d", model.name()),
            &short,
        ));
        nulls_long.push(threshold_rows(
            &mut report,
            &format!("{} noise, P 20–50 d", model.name()),
            &long,
        ));
    }
    let mut real: Vec<(String, Vec<Point>)> = fs::read_dir(&real_dir)
        .map(|entries| {
            entries
                .filter_map(|e| e.ok().map(|e| e.path()))
                .filter(|p| {
                    p.extension().is_some_and(|x| x == "csv")
                        && p.file_name().is_some_and(|f| f != "manifest.csv")
                })
                .map(|p| {
                    (
                        p.file_name().unwrap().to_string_lossy().into_owned(),
                        read_curve(&p),
                    )
                })
                .collect()
        })
        .unwrap_or_default();
    real.sort_by(|a, b| a.0.cmp(&b.0));
    real.truncate(real_limit);
    let real_cuts = (!real.is_empty()).then(|| {
        let inverted: Vec<[f64; 2]> = real
            .iter()
            .enumerate()
            .map(|(i, (name, points))| {
                let o = h.search(&invert(points), (1.0, 10.0), i < 2);
                record(&mut h, "null_real_inverted", "real", i, name, &o, &[]);
                o.top
            })
            .collect();
        let cuts = threshold_rows(
            &mut report,
            "inverted real (sector 32, 600 s), P 1–10 d",
            &inverted,
        );
        (cuts, inverted)
    });

    writeln!(
        report,
        "\n### Recovery at matched false-alarm rate (white → **correlated**)\n"
    )
    .unwrap();
    writeln!(report, "| Family | Injections | 10% FAR (top 12) | 5% FAR (top 12) | 1% FAR (top 12) | Top 1 @ 5% | Top 5 @ 5% | Top 12 @ 5% |").unwrap();
    writeln!(
        report,
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"
    )
    .unwrap();
    for (m, model) in MODELS.iter().enumerate() {
        // Varied durations at fixed white-noise SNR targets.
        for hours in [1.0, 2.0, 4.0, 6.0] {
            let mut trials = Vec::new();
            for snr in [7.0, 9.0, 12.0] {
                for rep in 0..n(12.0) {
                    let seed = 100_000 * (m as u64 + 1)
                        + (hours as u64) * 10_000
                        + (snr as u64) * 100
                        + rep as u64;
                    let mut rng = Rng(seed ^ 0xD1B5_4A32_D192_ED03);
                    let mut points = simulate(*model, seed);
                    let period = rng.range(3.0, 7.0);
                    let mut planet = Planet {
                        period,
                        epoch: rng.range(0.5, 0.5 + period),
                        duration: hours / 24.0,
                        depth: 0.0,
                    };
                    planet.depth = depth_for_snr(&points, &planet, snr, 0.0);
                    planet.inject(&mut points);
                    let o = h.search(&points, (1.0, 10.0), rep == 0 && snr == 7.0);
                    let span = o.reports[0].span_days;
                    let t = single(&o, |period, epoch, _| {
                        ephemeris_matches_planet(period, epoch, &planet, span)
                    });
                    record(
                        &mut h,
                        &format!("duration_{hours}h"),
                        model.name(),
                        rep,
                        &format!("snr{snr} P{period:.4} depth{:.1}", planet.depth * 1e6),
                        &o,
                        &[t],
                    );
                    trials.push(t);
                }
            }
            table_row(
                &mut report,
                &format!("{} noise, {hours} h transits (SNR 7/9/12)", model.name()),
                &trials,
                &nulls_short[m],
            );
        }
        // Two-event planets: one transit per sector, P from 20 to 50 days.
        let mut trials = Vec::new();
        for (k, snr) in [8.0, 10.0, 13.0].iter().enumerate() {
            let mut rep = 0;
            let mut seed = 300_000 * (m as u64 + 1) + k as u64 * 1000;
            while rep < n(15.0) {
                seed += 1;
                let mut rng = Rng(seed ^ 0x94D0_49BB_1331_11EB);
                let mut points = simulate(*model, seed);
                let (t1, t2) = (rng.range(0.5, 27.0), rng.range(90.5, 117.0));
                let cycles = ((t2 - t1) / rng.range(20.0, 50.0)).round().max(1.0);
                let period = (t2 - t1) / cycles;
                let mut planet = Planet {
                    period,
                    epoch: t1,
                    duration: 0.15,
                    depth: 0.0,
                };
                if !(20.0..=50.0).contains(&period) || planet.observed_events(&points) != 2 {
                    continue;
                }
                planet.depth = depth_for_snr(&points, &planet, *snr, 0.0);
                planet.inject(&mut points);
                let o = h.search(&points, (20.0, 50.0), rep == 0);
                let t = single(&o, |period, epoch, _| {
                    predicts_ephemeris(period, epoch, &[t1, t2], planet.duration)
                });
                record(
                    &mut h,
                    "two_event",
                    model.name(),
                    rep,
                    &format!("snr{snr} P{period:.3}"),
                    &o,
                    &[t],
                );
                trials.push(t);
                rep += 1;
            }
        }
        table_row(
            &mut report,
            &format!(
                "{} noise, two-event planets, 3.6 h (SNR 8/10/13)",
                model.name()
            ),
            &trials,
            &nulls_long[m],
        );
        // Two-planet systems: each planet is recovered if any listed dip matches it.
        let (mut inner, mut outer) = (Vec::new(), Vec::new());
        let mut seed = 400_000 * (m as u64 + 1);
        while inner.len() < n(40.0) {
            seed += 1;
            let mut rng = Rng(seed ^ 0xBF58_476D_1CE4_E5B9);
            let (p1, p2) = (rng.range(1.5, 3.5), rng.range(4.0, 9.5));
            if ((p2 / p1) - (p2 / p1).round()).abs() < 0.05 {
                continue;
            }
            let mut points = simulate(*model, seed);
            let mut a = Planet {
                period: p1,
                epoch: rng.range(0.5, 0.5 + p1),
                duration: 2.0 / 24.0,
                depth: 0.0,
            };
            let mut b = Planet {
                period: p2,
                epoch: rng.range(0.5, 0.5 + p2),
                duration: 3.0 / 24.0,
                depth: 0.0,
            };
            a.depth = depth_for_snr(&points, &a, 12.0, 0.0);
            b.depth = depth_for_snr(&points, &b, 10.0, 0.0);
            a.inject(&mut points);
            b.inject(&mut points);
            let o = h.search(&points, (1.0, 10.0), inner.is_empty());
            let span = o.reports[0].span_days;
            let any = |planet: &Planet| Trial {
                matched: [0, 1].map(|mode| {
                    dips(&o.reports[mode])
                        .iter()
                        .enumerate()
                        .find_map(|(index, c)| {
                            matched_ephemeris_score(c, planet, span).map(|score| Matched {
                                score,
                                rank: index + 1,
                            })
                        })
                }),
            };
            let (ta, tb) = (any(&a), any(&b));
            record(
                &mut h,
                "two_planet",
                model.name(),
                inner.len(),
                &format!("P{p1:.3}/{p2:.3}"),
                &o,
                &[ta, tb],
            );
            inner.push(ta);
            outer.push(tb);
        }
        table_row(
            &mut report,
            &format!("{} noise, two-planet: inner (2 h, SNR 12)", model.name()),
            &inner,
            &nulls_short[m],
        );
        table_row(
            &mut report,
            &format!("{} noise, two-planet: outer (3 h, SNR 10)", model.name()),
            &outer,
            &nulls_short[m],
        );
    }
    if let Some((cuts, star_nulls)) = &real_cuts {
        // Variable stars dominate pooled real nulls, and their own signals take
        // the top rank, so a real injection counts when any listed dip matches.
        // Per-star matching compares it with the same star's inverted curve.
        let mut per_star_lines = Vec::new();
        for hours in [2.0, 4.0] {
            let mut trials = Vec::new();
            let mut per_star = [0, 0];
            for (i, (name, curve)) in real.iter().enumerate() {
                for snr in [8.0, 12.0] {
                    let seed = 700_000 + (hours as u64) * 1000 + i as u64 * 10 + snr as u64;
                    let mut rng = Rng(seed);
                    let mut points = curve.clone();
                    let period = rng.range(1.5, 8.0);
                    let mut planet = Planet {
                        period,
                        epoch: points[0].time + rng.range(0.2, period),
                        duration: hours / 24.0,
                        depth: 0.0,
                    };
                    planet.depth = depth_for_snr(&points, &planet, snr, sample_sigma(curve));
                    planet.inject(&mut points);
                    let o = h.search(&points, (1.0, 10.0), false);
                    let span = o.reports[0].span_days;
                    let t = Trial {
                        matched: [0, 1].map(|mode| {
                            dips(&o.reports[mode])
                                .iter()
                                .enumerate()
                                .find_map(|(index, c)| {
                                    matched_ephemeris_score(c, &planet, span).map(|score| Matched {
                                        score,
                                        rank: index + 1,
                                    })
                                })
                        }),
                    };
                    for (mode, count) in per_star.iter_mut().enumerate() {
                        if t.matched[mode].is_some_and(|m| m.score > star_nulls[i][mode]) {
                            *count += 1;
                        }
                    }
                    record(
                        &mut h,
                        &format!("real_{hours}h"),
                        "real",
                        i,
                        &format!("{name} snr{snr} P{period:.4}"),
                        &o,
                        &[t],
                    );
                    trials.push(t);
                }
            }
            table_row(
                &mut report,
                &format!("real curves, {hours} h transits (SNR 8/12), pooled nulls"),
                &trials,
                cuts,
            );
            per_star_lines.push(format!(
                "- Real curves, {hours} h: {} of {} injections beat the same star's inverted-curve top score with the white statistic, **{}** with the correlated one.",
                per_star[0],
                trials.len(),
                per_star[1]
            ));
        }
        writeln!(report, "\n{}", per_star_lines.join("\n")).unwrap();
    }
    writeln!(
        report,
        "\nSearch time: white {:.0} s, correlated {:.0} s. Serial/parallel checks: {} curves, {} mismatches. Total {:.0} s.",
        h.seconds[0], h.seconds[1], h.determinism_checks, h.determinism_failures, started.elapsed().as_secs_f64()
    )
    .unwrap();
    fs::create_dir_all(&out_dir).expect("output directory");
    fs::write(out_dir.join("raw.csv"), &h.raw).expect("raw results");
    fs::write(out_dir.join("summary.md"), &report).expect("summary");
    print!("{report}");
}

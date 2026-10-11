use astrosearch::{
    analyze_file, provenance_json, render_csv, render_csv_with_tier, CandidateGrouping, Config,
    TrendMethod,
};
use std::{
    env, fs,
    io::{self, Write},
    path::{Path, PathBuf},
    process,
    sync::atomic::{AtomicU64, Ordering},
};

static TEMP_SEQUENCE: AtomicU64 = AtomicU64::new(0);

fn main() {
    if let Err(error) = run() {
        eprintln!("astrosearch: {error}");
        process::exit(2);
    }
}

fn run() -> Result<(), String> {
    let raw_args: Vec<String> = env::args().skip(1).collect();
    let has_fixed_period_grid = raw_args.iter().any(|arg| arg == "--periods");
    let mut selected_profile = String::from("fast");
    let mut index = 0;
    while index < raw_args.len() {
        if raw_args[index] == "--search-profile" {
            let Some(profile) = raw_args.get(index + 1) else {
                return Err("missing value after --search-profile".into());
            };
            if !matches!(profile.as_str(), "fast" | "standard") {
                return Err(format!("unknown search profile {profile}; use fast or standard"));
            }
            selected_profile = profile.clone();
            index += 2;
        } else {
            index += 1;
        }
    }
    let mut args = raw_args.into_iter();
    let mut input = None;
    let mut output = None;
    let mut periodogram = None;
    let mut config = if selected_profile == "standard" {
        Config::standard()
    } else {
        Config::profile_j()
    };
    while let Some(arg) = args.next() {
        let value = |args: &mut std::vec::IntoIter<String>| {
            args.next()
                .ok_or_else(|| format!("missing value after {arg}"))
        };
        match arg.as_str() {
            "--search-profile" => {
                let _ = value(&mut args)?;
            }
            "--input" | "-i" => input = Some(value(&mut args)?),
            "--output" | "-o" => output = Some(value(&mut args)?),
            "--min-period" => config.min_period = parse_num(&value(&mut args)?, "min-period")?,
            "--max-period" => config.max_period = parse_num(&value(&mut args)?, "max-period")?,
            "--stellar-density" => {
                config.stellar_density_g_cm3 =
                    Some(parse_num(&value(&mut args)?, "stellar-density")?)
            }
            "--physical-duration-floor" => config.physical_duration_floor = true,
            "--two-stage" => config.two_stage_factor = parse_num(&value(&mut args)?, "two-stage")?,
            "--coarse-ffa" => config.coarse_ffa = true,
            "--stop-passes-below-detection" => config.stop_passes_below_detection = true,
            "--two-stage-keep" => {
                config.two_stage_keep = parse_num(&value(&mut args)?, "two-stage-keep")?
            }
            "--bin-minutes" => {
                config.bin_days = parse_num::<f64>(&value(&mut args)?, "bin-minutes")? / 1440.0
            }
            "--periods" => config.periods = parse_num(&value(&mut args)?, "periods")?,
            "--durations" => config.durations = parse_num(&value(&mut args)?, "durations")?,
            "--trend-days" => config.trend_days = parse_num(&value(&mut args)?, "trend-days")?,
            "--segment-gap-days" => {
                config.segment_gap_days = parse_num(&value(&mut args)?, "segment-gap-days")?
            }
            "--trend-method" => config.trend_method = value(&mut args)?.parse::<TrendMethod>()?,
            "--transit-search-passes" => {
                config.transit_search_passes =
                    parse_num(&value(&mut args)?, "transit-search-passes")?
            }
            "--invert" => config.invert = true,
            "--min-event-points" => {
                config.min_event_points = parse_num(&value(&mut args)?, "min-event-points")?
            }
            "--threads" => config.threads = parse_num(&value(&mut args)?, "threads")?,
            "--scoring" => config.scoring = value(&mut args)?.parse()?,
            "--candidate-grouping" => {
                config.candidate_grouping = value(&mut args)?.parse::<CandidateGrouping>()?
            }
            "--min-transit-score" => {
                config.min_transit_score = parse_num(&value(&mut args)?, "min-transit-score")?
            }
            "--max-period-trials" => {
                config.max_period_trials = parse_num(&value(&mut args)?, "max-period-trials")?
            }
            "--detection-threshold" => {
                config.detection_threshold = parse_num(&value(&mut args)?, "detection-threshold")?
            }
            "--transit-candidates" => {
                config.transit_candidates = parse_num(&value(&mut args)?, "transit-candidates")?
            }
            "--periodogram" => {
                periodogram = Some(value(&mut args)?);
                config.keep_periodogram = true;
            }
            "--help" | "-h" => {
                help();
                return Ok(());
            }
            other => return Err(format!("unknown option {other}; use --help")),
        }
    }
    if selected_profile == "fast" && has_fixed_period_grid {
        return Err("--periods fixes the grid; use --search-profile standard with --periods".into());
    }
    let input = input.ok_or("--input is required")?;
    for path in output.iter().chain(&periodogram) {
        let canonical_input = fs::canonicalize(&input).map_err(|e| format!("open {input}: {e}"))?;
        if fs::canonicalize(path).is_ok_and(|p| p == canonical_input) || same_file(&input, path) {
            return Err(
                "output must differ from input; refusing to overwrite the source light curve"
                    .into(),
            );
        }
    }
    if let (Some(candidates), Some(trials)) = (&output, &periodogram) {
        if candidates == trials
            || same_destination(candidates, trials)
            || same_file(candidates, trials)
        {
            return Err("candidate output and periodogram output must be different files".into());
        }
    }
    let report = analyze_file(&input, config)?;
    if let Some(path) = &periodogram {
        // One row per trial period with an accepted window, in period order.
        atomic_write(path, render_csv(&report.periodogram).as_bytes())?;
    }
    let csv = render_csv_with_tier(&report.candidates, config.detection_threshold);
    if let Some(path) = output {
        atomic_write(&path, csv.as_bytes())?;
        // Version, configuration, and input digest that produced this file.
        let bytes = fs::read(&input).map_err(|e| format!("read {input}: {e}"))?;
        let record = provenance_json(&input, &bytes, &config, &report);
        atomic_write(&format!("{path}.provenance.json"), record.as_bytes())?;
    } else {
        let result = io::stdout().lock().write_all(csv.as_bytes());
        if let Err(error) = result {
            if error.kind() != io::ErrorKind::BrokenPipe {
                return Err(format!("write stdout: {error}"));
            }
        }
    }
    for warning in &report.warnings {
        eprintln!("warning: {warning}");
    }
    eprintln!("transit period trials: {}", report.transit_trials);
    for (i, s) in report.segments.iter().enumerate() {
        let beta: Vec<String> = s
            .red_noise
            .iter()
            .map(|(days, beta)| format!("{:.1}h:{beta:.2}", days * 24.0))
            .collect();
        eprintln!(
            "segment {}: {:.4}-{:.4}, {} observations, cadence {:.1} s, noise {:.0} ppm, beta {}",
            i + 1,
            s.start_time,
            s.end_time,
            s.observations,
            s.cadence_days * 86400.0,
            s.noise_ppm,
            beta.join(" ")
        );
    }
    eprintln!(
        "searched {} valid observations across {:.3} days; {} ranked candidates",
        report.observations,
        report.span_days,
        report.candidates.len()
    );
    Ok(())
}

fn atomic_write(path: &str, contents: &[u8]) -> Result<(), String> {
    let destination = Path::new(path);
    let parent = destination
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let name = destination
        .file_name()
        .ok_or_else(|| format!("invalid output path: {path}"))?;
    let mut temp = PathBuf::new();
    let mut file = None;
    for _ in 0..100 {
        let sequence = TEMP_SEQUENCE.fetch_add(1, Ordering::Relaxed);
        temp = parent.join(format!(
            ".{}.{}.{}.tmp",
            name.to_string_lossy(),
            process::id(),
            sequence
        ));
        match fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&temp)
        {
            Ok(created) => {
                file = Some(created);
                break;
            }
            Err(error) if error.kind() == io::ErrorKind::AlreadyExists => continue,
            Err(error) => return Err(format!("create temporary output for {path}: {error}")),
        }
    }
    let mut file = file.ok_or_else(|| format!("could not allocate temporary output for {path}"))?;
    let result = (|| {
        file.write_all(contents)
            .map_err(|e| format!("write temporary output for {path}: {e}"))?;
        file.sync_all()
            .map_err(|e| format!("sync temporary output for {path}: {e}"))?;
        drop(file);
        fs::rename(&temp, destination).map_err(|e| format!("publish output {path}: {e}"))
    })();
    if result.is_err() {
        let _ = fs::remove_file(&temp);
    }
    result
}

fn same_destination(a: &str, b: &str) -> bool {
    fn canonical_destination(path: &str) -> Option<PathBuf> {
        let path = Path::new(path);
        let name = path.file_name()?;
        let parent = path
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or(Path::new("."));
        Some(fs::canonicalize(parent).ok()?.join(name))
    }
    canonical_destination(a)
        .zip(canonical_destination(b))
        .is_some_and(|(a, b)| a == b)
}

#[cfg(unix)]
fn same_file(input: &str, output: &str) -> bool {
    use std::os::unix::fs::MetadataExt;
    match (fs::metadata(input), fs::metadata(output)) {
        (Ok(a), Ok(b)) => a.dev() == b.dev() && a.ino() == b.ino(),
        _ => false,
    }
}

#[cfg(not(unix))]
fn same_file(_: &str, _: &str) -> bool {
    false
}

fn parse_num<T: std::str::FromStr>(value: &str, name: &str) -> Result<T, String> {
    value
        .parse()
        .map_err(|_| format!("invalid {name}: {value}"))
}

fn help() {
    println!("astrosearch — TESS-style light-curve candidate search\n\nUSAGE:\n  astrosearch --input lightcurve.csv [--output candidates.csv]\n              [--search-profile fast|standard] [--min-period DAYS] [--max-period DAYS]\n              [--periods N] [--stellar-density G_PER_CM3]\n              [--durations N] [--trend-days DAYS] [--segment-gap-days DAYS]\n              [--trend-method median|biweight|none] [--transit-search-passes N]\n              [--min-event-points N] [--threads N] [--scoring correlated|white]\n              [--candidate-grouping ephemeris|events] [--min-transit-score X]\n              [--detection-threshold X] [--max-period-trials N]\n              [--bin-minutes M] [--physical-duration-floor]\n              [--two-stage F] [--two-stage-keep K] [--coarse-ffa]\n              [--stop-passes-below-detection]\n              [--transit-candidates N] [--periodogram trials.csv] [--invert]\n\nInput needs a time/BTJD/BJD and positive-flux column; errors and quality are optional.\nCSV or TSV, including quoted fields, supported. Period bounds apply to both periodic detectors.\nProfile J is the default: 10-minute binning, physical duration floor, 7x coarse FFA\nsearch, 2-million period-trial cap, five passes, and QLP-calibrated jackknife tier\n25.36783. Apply that tier only to the matching QLP search contract; use a matching\nproduct-specific calibration for other data. `--search-profile standard` selects the\nfull-grid baseline and legacy tier; `--periods` requires this profile.\nThreads default to 0 (automatic); use 1 for serial execution.\nScoring defaults to correlated (per-segment, duration-dependent noise); white is the\noriginal independent-error statistic. Both appear in the output (score, white_score).\n--stellar-density uses g/cm^3 and assumes b <= 0.95 and e <= 0.9 for a target-specific\nduration floor; omit it for the cadence-only completeness search.\nTrend defaults to a sliding-window Tukey biweight with a 1.05-day window (three times the 0.35-day maximum searched duration).\n--candidate-grouping defaults to events: aliases sharing observed transit events are grouped and kept in the CSV;\nephemeris is the earlier period-drift deduplication.\n--transit-search-passes (default 5) masks each pass's strongest event windows and reruns the transit search;\nruntime grows with passes, and 1 disables masking.\n--bin-minutes M overrides Profile J's 10-minute binning; 0 disables binning.\n--physical-duration-floor uses a geometry-informed minimum duration; this saves trials,\nbut strongly grazing or eccentric short transits may not be searched.\n--two-stage F overrides Profile J's factor of 7; --two-stage-keep defaults to 512.\n--coarse-ffa ranks the cells with the fast folding algorithm over all periods.\n--stop-passes-below-detection skips further passes once a pass finds nothing at the detection tier.\n--max-period-trials defaults to 2000000 for Profile J and 250000 for standard.\nLong baselines may still reach the cap; warnings report a shortfall.\n--min-transit-score (default 6) is only a screening floor; many inverted null controls exceed it.\njackknife_score is the score without the strongest single event; the tier column labels\ntransit candidates above the configured threshold as detection, others as screening.\n--invert reflects flux about its median for an inverted-flux systematic control.\nEpochs are window centers; scores rank unvalidated light-curve triggers.");
}

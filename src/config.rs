//! Search options, their defaults, and validation.

/// Default cap on trial periods per search for Profile J.
pub const MAX_PERIOD_TRIALS: usize = 2_000_000;
/// Grid cap used by the explicit standard profile.
pub const STANDARD_MAX_PERIOD_TRIALS: usize = 250_000;
/// Hard ceiling on `Config::max_period_trials` (memory and runtime guard).
pub const PERIOD_TRIALS_CEILING: usize = 20_000_000;
/// Calibrated five-pass Profile J detection tier.
pub const DETECTION_THRESHOLD: f64 = 25.36783;
/// Detection tier for the explicit standard profile's legacy calibration.
pub const STANDARD_DETECTION_THRESHOLD: f64 = 23.302528953788222;

/// Statistic that selects and ranks periodic-dip windows.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Scoring {
    /// Box SNR assuming independent errors (the original statistic).
    White,
    /// Box SNR whose in-window variance is inflated by each segment's
    /// correlated-noise factor at the window's duration.
    Correlated,
}

impl std::str::FromStr for Scoring {
    type Err = String;
    fn from_str(value: &str) -> Result<Self, String> {
        match value {
            "white" => Ok(Self::White),
            "correlated" => Ok(Self::Correlated),
            other => Err(format!("unknown scoring {other}; use white or correlated")),
        }
    }
}

/// Rule used to group periodic-dip candidates before applying the output cap.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CandidateGrouping {
    /// Group only candidates with nearly identical ephemerides (legacy behavior).
    Ephemeris,
    /// Group candidates that explain substantially the same observed events.
    Events,
}

impl std::str::FromStr for CandidateGrouping {
    type Err = String;
    fn from_str(value: &str) -> Result<Self, String> {
        match value {
            "ephemeris" => Ok(Self::Ephemeris),
            "events" => Ok(Self::Events),
            other => Err(format!(
                "unknown candidate-grouping {other}; use ephemeris or events"
            )),
        }
    }
}

/// Robust trend estimator used independently in each gap/cadence segment.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TrendMethod {
    Median,
    Biweight,
    /// No detrending beyond one median level per segment, for light curves
    /// already flattened elsewhere (for example by a Gaussian-process fit).
    None,
}

impl std::str::FromStr for TrendMethod {
    type Err = String;
    fn from_str(value: &str) -> Result<Self, String> {
        match value {
            "median" => Ok(Self::Median),
            "biweight" => Ok(Self::Biweight),
            "none" => Ok(Self::None),
            other => Err(format!(
                "unknown trend-method {other}; use median, biweight, or none"
            )),
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Config {
    pub min_period: f64,
    pub max_period: f64,
    /// Optional stellar mean density in g/cm^3. When supplied, the transit
    /// search skips durations shorter than a conservative geometry bound.
    pub stellar_density_g_cm3: Option<f64>,
    /// Use half the circular central-transit duration (stellar density, or
    /// solar) as the shortest searched duration instead of the cadence floor.
    /// Needs far fewer trial periods on long baselines; very short (grazing or
    /// eccentric) transits are then not searched.
    pub physical_duration_floor: bool,
    /// Bin samples to this width (days) before detrending; 0 disables. Bins
    /// stay within gaps and slower cadences pass through unchanged.
    pub bin_days: f64,
    /// Two-stage search: evaluate one period per this many fine grid periods,
    /// then refine the best cells at full resolution (0 disables).
    pub two_stage_factor: usize,
    /// Coarse cells refined in a two-stage search.
    pub two_stage_keep: usize,
    /// Skip further masked passes once a pass's best transit score is below
    /// `detection_threshold`: later passes could only add screening-tier
    /// candidates.
    pub stop_passes_below_detection: bool,
    /// Rank coarse cells with the fast folding algorithm (all periods at
    /// duty-cycle-scaled resolution) instead of one grid period per cell.
    pub coarse_ffa: bool,
    /// Minimum number of period/frequency trials; resolution can require more.
    pub periods: usize,
    pub durations: usize,
    pub trend_days: f64,
    /// Split segments across gaps longer than this duration (days).
    pub segment_gap_days: f64,
    pub trend_method: TrendMethod,
    pub min_event_points: usize,
    /// Worker count; zero automatically uses available CPU parallelism.
    pub threads: usize,
    pub scoring: Scoring,
    /// Event-based alias grouping is selectable; ephemeris grouping is the baseline.
    pub candidate_grouping: CandidateGrouping,
    /// Periodic-dip windows below this score are discarded before vetting.
    /// This is a screening floor: 97% of inverted null controls in the
    /// known-answer benchmark produce a candidate above it.
    pub min_transit_score: f64,
    /// Transit candidates whose jackknife score (`Candidate::jackknife_score`)
    /// reaches this value are labeled `detection`, the rest `screening`. The
    /// Profile J default is the deployment tier from the verified five-pass
    /// QLP calibration; product and search contract must match before treating
    /// this label as calibrated.
    pub detection_threshold: f64,
    /// Cap on trial periods. When the resolution rule needs more, the grid is
    /// coarsened and a warning is emitted; raising it costs runtime roughly
    /// in proportion (see docs/GRID_CAP_STUDY.md).
    pub max_period_trials: usize,
    /// Distinct periodic-dip ephemerides kept after deduplication.
    pub transit_candidates: usize,
    /// Maximum number of sequential masked transit searches.
    pub transit_search_passes: usize,
    /// Reflect flux around its median for an inverted-curve false-alarm control.
    pub invert: bool,
    /// Return each trial period's best accepted window in `Report::periodogram`.
    pub keep_periodogram: bool,
}

impl Default for Config {
    fn default() -> Self {
        Self::profile_j()
    }
}

impl Config {
    /// Recommended general-purpose configuration: five-pass Profile J.
    /// The threshold is the deployment tier from the verified QLP calibration
    /// `profile-j-fivepass-20261010-fa05-verified.json`.
    pub fn profile_j() -> Self {
        Self {
            min_period: 0.15,
            max_period: 30.0,
            stellar_density_g_cm3: None,
            physical_duration_floor: true,
            bin_days: 10.0 / 1440.0,
            two_stage_factor: 7,
            two_stage_keep: 512,
            coarse_ffa: true,
            stop_passes_below_detection: false,
            periods: 2400,
            durations: 10,
            trend_days: 1.05,
            segment_gap_days: 0.5,
            trend_method: TrendMethod::Biweight,
            min_event_points: 2,
            threads: 0,
            scoring: Scoring::Correlated,
            candidate_grouping: CandidateGrouping::Events,
            min_transit_score: 6.0,
            detection_threshold: DETECTION_THRESHOLD,
            max_period_trials: MAX_PERIOD_TRIALS,
            transit_candidates: 12,
            transit_search_passes: 5,
            invert: false,
            keep_periodogram: false,
        }
    }

    /// Original full-grid configuration for explicit baseline comparisons.
    pub fn standard() -> Self {
        Self {
            physical_duration_floor: false,
            bin_days: 0.0,
            two_stage_factor: 0,
            coarse_ffa: false,
            stop_passes_below_detection: false,
            detection_threshold: STANDARD_DETECTION_THRESHOLD,
            max_period_trials: STANDARD_MAX_PERIOD_TRIALS,
            transit_search_passes: 1,
            ..Self::profile_j()
        }
    }

    pub fn validate(&self) -> Result<(), String> {
        if ![
            self.min_period,
            self.max_period,
            self.trend_days,
            self.segment_gap_days,
        ]
        .iter()
        .all(|v| v.is_finite() && *v > 0.0 && (1.0 / v).is_finite())
            || self.max_period <= self.min_period
        {
            return Err("period bounds and trend-days must be finite, positive numbers; max-period must exceed min-period".into());
        }
        if self
            .stellar_density_g_cm3
            .is_some_and(|rho| !rho.is_finite() || rho <= 0.0)
        {
            return Err("stellar-density must be a finite, positive value in g/cm^3".into());
        }
        if self.two_stage_factor == 1 || self.two_stage_factor > 1024 {
            return Err("two-stage factor must be 0 (off) or between 2 and 1024".into());
        }
        if self.coarse_ffa && self.two_stage_factor < 2 {
            return Err("coarse-ffa needs --two-stage".into());
        }
        if !(1..=1_000_000).contains(&self.two_stage_keep) {
            return Err("two-stage-keep must be between 1 and 1000000".into());
        }
        if !(self.bin_days.is_finite() && (0.0..=1.0).contains(&self.bin_days)) {
            return Err("bin-minutes must be between 0 and 1440".into());
        }
        if !(10..=PERIOD_TRIALS_CEILING).contains(&self.max_period_trials) {
            return Err(format!(
                "max-period-trials must be between 10 and {PERIOD_TRIALS_CEILING}"
            ));
        }
        if !(10..=self.max_period_trials).contains(&self.periods) {
            return Err(format!(
                "periods must be between 10 and max-period-trials ({})",
                self.max_period_trials
            ));
        }
        if !(1..=128).contains(&self.durations) {
            return Err("durations must be between 1 and 128".into());
        }
        if !(2..=1024).contains(&self.min_event_points) {
            return Err("min-event-points must be between 2 and 1024".into());
        }
        if self.threads > 256 {
            return Err("threads must be between 0 (automatic) and 256".into());
        }
        if !(self.min_transit_score.is_finite() && self.min_transit_score > 0.0) {
            return Err("min-transit-score must be a finite, positive number".into());
        }
        if !(self.detection_threshold.is_finite() && self.detection_threshold > 0.0) {
            return Err("detection-threshold must be a finite, positive number".into());
        }
        if !(1..=100_000).contains(&self.transit_candidates) {
            return Err("transit-candidates must be between 1 and 100000".into());
        }
        if !(1..=32).contains(&self.transit_search_passes) {
            return Err("transit-search-passes must be between 1 and 32".into());
        }
        Ok(())
    }
}

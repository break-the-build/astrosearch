//! Search options, their defaults, and validation.

pub const MAX_PERIOD_TRIALS: usize = 250_000;

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
    pub min_transit_score: f64,
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
        Self {
            min_period: 0.15,
            max_period: 30.0,
            stellar_density_g_cm3: None,
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
            transit_candidates: 12,
            transit_search_passes: 5,
            invert: false,
            keep_periodogram: false,
        }
    }
}

impl Config {
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
        if !(10..=MAX_PERIOD_TRIALS).contains(&self.periods) {
            return Err(format!(
                "periods must be between 10 and {MAX_PERIOD_TRIALS}"
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
        if !(1..=100_000).contains(&self.transit_candidates) {
            return Err("transit-candidates must be between 1 and 100000".into());
        }
        if !(1..=32).contains(&self.transit_search_passes) {
            return Err("transit-search-passes must be between 1 and 32".into());
        }
        Ok(())
    }
}

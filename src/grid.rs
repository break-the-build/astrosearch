//! Trial-period grid, trial budgets, and the minimum searched transit duration.

/// Requested trial count (`Config::periods`) and its cap
/// (`Config::max_period_trials`).
#[derive(Clone, Copy)]
pub(crate) struct TrialBudget {
    pub(crate) requested: usize,
    pub(crate) cap: usize,
}

impl TrialBudget {
    pub(crate) fn from_config(config: &crate::Config) -> Self {
        Self {
            requested: config.periods,
            cap: config.max_period_trials,
        }
    }
}

/// Trials needed for the target resolution, at least `budget.requested`,
/// capped at `budget.cap` with a warning when the cap loses resolution.
pub(crate) fn trial_count(
    required: usize,
    budget: TrialBudget,
    detector: &str,
    warnings: &mut Vec<String>,
) -> usize {
    let cap = budget.cap;
    if required > cap {
        warnings.push(format!("{detector} search needs {required} trials for target resolution; capped at {cap}, so narrow peaks may be missed"));
    }
    budget.requested.max(required).clamp(2, cap)
}

/// Solar mean density, the physical floor's default when none is supplied.
pub(crate) const SOLAR_DENSITY_G_CM3: f64 = 1.41;
/// Physical floor as a fraction of the central-transit duration (as in the
/// common "grazing factor" prior; b = 0.87 for a small planet).
pub(crate) const PHYSICAL_FLOOR_FRACTION: f64 = 0.5;

/// How the shortest searched transit duration is set for each trial period.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct DurationModel {
    pub(crate) density_g_cm3: Option<f64>,
    pub(crate) physical: bool,
}

impl DurationModel {
    pub(crate) fn from_config(config: &crate::Config) -> Self {
        Self {
            density_g_cm3: config.stellar_density_g_cm3,
            physical: config.physical_duration_floor,
        }
    }

    #[cfg(test)]
    pub(crate) const CADENCE: Self = Self {
        density_g_cm3: None,
        physical: false,
    };

    fn cadence_only(self) -> bool {
        self.density_g_cm3.is_none() && !self.physical
    }
}

fn scaled_axis(period: f64, density_g_cm3: f64) -> f64 {
    const G: f64 = 6.67430e-11;
    let period_seconds = period * 86_400.0;
    (G * density_g_cm3 * 1_000.0 * period_seconds.powi(2) / (3.0 * std::f64::consts::PI)).cbrt()
}

/// Shortest searched duration at `period`:
/// * cadence floor (2 x cadence) by default;
/// * with a stellar density, a conservative floor assuming b <= 0.95 and
///   e <= 0.9 (explicit search priors, not universal limits);
/// * with `physical`, PHYSICAL_FLOOR_FRACTION of the circular central-transit
///   duration (solar density unless given). Much shorter transits, such as
///   strongly grazing or highly eccentric ones, are then not searched.
pub(crate) fn minimum_duration(period: f64, cadence: f64, model: DurationModel) -> f64 {
    let sampled = 2.0 * cadence;
    if model.physical {
        let density = model.density_g_cm3.unwrap_or(SOLAR_DENSITY_G_CM3);
        let central = period / std::f64::consts::PI
            * (1.0 / scaled_axis(period, density)).clamp(0.0, 1.0).asin();
        return sampled.max(PHYSICAL_FLOOR_FRACTION * central);
    }
    let Some(density) = model.density_g_cm3 else {
        return sampled;
    };
    const MAX_IMPACT: f64 = 0.95;
    const MAX_ECCENTRICITY: f64 = 0.9;
    let chord = (1.0_f64 - MAX_IMPACT.powi(2)).sqrt();
    let circular = period / std::f64::consts::PI
        * (chord / scaled_axis(period, density))
            .clamp(0.0, 1.0)
            .asin();
    let eccentric_shortening = (1.0 - MAX_ECCENTRICITY.powi(2)).sqrt() / (1.0 + MAX_ECCENTRICITY);
    sampled.max(circular * eccentric_shortening)
}

/// Build a target-specific, duty-cycle-aware log-period grid. Adjacent samples
/// are spaced so accumulated timing drift stays within half the shortest
/// searched duration. A left-Riemann integral over log-period conservatively
/// estimates the required grid because the physical duration floor is
/// nondecreasing with period.
pub(crate) fn adaptive_period_grid(
    min_period: f64,
    max_period: f64,
    span: f64,
    cadence: f64,
    model: DurationModel,
    budget: TrialBudget,
    warnings: &mut Vec<String>,
) -> Vec<f64> {
    const INTEGRATION_STEPS: usize = 4096;
    let log_min = min_period.ln();
    let log_range = max_period.ln() - log_min;
    if model.cadence_only() {
        let duration = minimum_duration(min_period, cadence, model);
        let max_log_step = (duration / (2.0 * span)).ln_1p();
        let required = ((log_range / max_log_step).ceil() as usize).saturating_add(1);
        let trials = trial_count(required, budget, "transit", warnings);
        let step = log_range / (trials - 1) as f64;
        return (0..trials)
            .map(|i| (log_min + step * i as f64).exp())
            .collect();
    }
    let log_step = log_range / INTEGRATION_STEPS as f64;
    let margin = 1.02;
    let mut resolution = Vec::with_capacity(INTEGRATION_STEPS + 1);
    resolution.push(0.0);
    for i in 0..INTEGRATION_STEPS {
        let period = (log_min + i as f64 * log_step).exp();
        let duration = minimum_duration(period, cadence, model);
        let max_log_step = (duration / (2.0 * span)).ln_1p();
        resolution.push(resolution[i] + margin * log_step / max_log_step);
    }

    // The physical duration curve uses a small margin for integration error.
    let total_resolution = resolution[INTEGRATION_STEPS];
    let required = (total_resolution.ceil() as usize).saturating_add(1);
    let trials = trial_count(required, budget, "transit", warnings);
    let mut grid = Vec::with_capacity(trials);
    for i in 0..trials {
        let target = total_resolution * i as f64 / (trials - 1) as f64;
        let upper = resolution
            .partition_point(|&value| value < target)
            .min(INTEGRATION_STEPS);
        let lower = upper.saturating_sub(1);
        let q0 = resolution[lower];
        let q1 = resolution[upper];
        let fraction = if q1 > q0 {
            ((target - q0) / (q1 - q0)).clamp(0.0, 1.0)
        } else {
            0.0
        };
        let log_period = log_min + (lower as f64 + fraction) * log_step;
        grid.push(log_period.exp().clamp(min_period, max_period));
    }
    if let Some(first) = grid.first_mut() {
        *first = min_period;
    }
    if let Some(last) = grid.last_mut() {
        *last = max_period;
    }
    grid
}

#[cfg(test)]
mod tests {
    use super::*;

    const TEST_BUDGET: TrialBudget = TrialBudget {
        requested: 10,
        cap: 250_000,
    };

    #[test]
    fn physical_floor_is_half_the_central_duration_and_shrinks_the_grid() {
        let physical = DurationModel {
            density_g_cm3: None,
            physical: true,
        };
        // Sun-like star, 365.25 d: central duration ~13 h; floor ~6.5 h.
        let floor = minimum_duration(365.25, 2.0 / 1440.0, physical);
        assert!((floor * 24.0 - 6.5).abs() < 0.3, "{floor}");
        // Never below the sampling floor.
        assert!(minimum_duration(0.2, 0.05, physical) >= 0.1);
        let budget = TrialBudget {
            requested: 10,
            cap: 20_000_000,
        };
        let span = 1800.0;
        let cadence = 600.0 / 86400.0;
        let cadence_grid = adaptive_period_grid(
            0.15,
            30.0,
            span,
            cadence,
            DurationModel::CADENCE,
            budget,
            &mut Vec::new(),
        );
        let physical_grid =
            adaptive_period_grid(0.15, 30.0, span, cadence, physical, budget, &mut Vec::new());
        assert!(
            physical_grid.len() * 3 < cadence_grid.len(),
            "{} vs {}",
            physical_grid.len(),
            cadence_grid.len()
        );
        // The drift rule holds for the physical floor at every step.
        for pair in physical_grid.windows(2) {
            let step = (pair[1] / pair[0]).ln();
            let duration = minimum_duration(pair[0], cadence, physical);
            assert!(step <= 1.03 * (duration / (2.0 * span)).ln_1p());
        }
    }

    #[test]
    fn adaptive_grid_preserves_log_spacing_without_density_priors() {
        let grid = adaptive_period_grid(
            0.5,
            20.0,
            1000.0,
            2.0 / 1440.0,
            DurationModel::CADENCE,
            TEST_BUDGET,
            &mut Vec::new(),
        );
        assert!(grid.len() > 10);
        assert_eq!(grid.first().copied(), Some(0.5));
        assert!((grid.last().copied().unwrap() - 20.0).abs() < 1e-12);
        let steps: Vec<f64> = grid
            .windows(2)
            .map(|pair| (pair[1] / pair[0]).ln())
            .collect();
        let mean = steps.iter().sum::<f64>() / steps.len() as f64;
        assert!(steps.iter().all(|step| (step / mean - 1.0).abs() < 1e-8));
    }

    #[test]
    fn density_grid_coarsens_with_period_and_respects_drift_bound() {
        let span = 10.0;
        let cadence = 2.0 / 86400.0;
        let grid = adaptive_period_grid(
            0.5,
            100.0,
            span,
            cadence,
            DurationModel {
                density_g_cm3: Some(1.0),
                physical: false,
            },
            TEST_BUDGET,
            &mut Vec::new(),
        );
        assert!(grid.len() > 10);
        let first_log_step = (grid[1] / grid[0]).ln();
        let last = grid.len() - 1;
        let last_log_step = (grid[last] / grid[last - 1]).ln();
        assert!(last_log_step > first_log_step);
        for pair in grid.windows(2) {
            let duration = minimum_duration(
                pair[0],
                cadence,
                DurationModel {
                    density_g_cm3: Some(1.0),
                    physical: false,
                },
            );
            let timing_drift = span * (pair[1] / pair[0] - 1.0);
            assert!(
                timing_drift <= duration * 0.51,
                "{pair:?}: {timing_drift} > {duration}"
            );
        }
    }
}

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

/// Sampling floor, or a stellar-density duration floor when target metadata
/// is supplied. The latter assumes b <= 0.95 and e <= 0.9; those are explicit
/// search priors, not universal limits on real transits.
pub(crate) fn minimum_duration(period: f64, cadence: f64, density_g_cm3: Option<f64>) -> f64 {
    let sampled = 2.0 * cadence;
    let Some(density) = density_g_cm3 else {
        return sampled;
    };
    const G: f64 = 6.67430e-11;
    const SECONDS_PER_DAY: f64 = 86_400.0;
    const MAX_IMPACT: f64 = 0.95;
    const MAX_ECCENTRICITY: f64 = 0.9;
    let period_seconds = period * SECONDS_PER_DAY;
    let density_kg_m3 = density * 1_000.0;
    let scaled_axis =
        (G * density_kg_m3 * period_seconds.powi(2) / (3.0 * std::f64::consts::PI)).cbrt();
    let chord = (1.0_f64 - MAX_IMPACT.powi(2)).sqrt();
    let circular = period / std::f64::consts::PI * (chord / scaled_axis).clamp(0.0, 1.0).asin();
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
    density_g_cm3: Option<f64>,
    budget: TrialBudget,
    warnings: &mut Vec<String>,
) -> Vec<f64> {
    const INTEGRATION_STEPS: usize = 4096;
    let log_min = min_period.ln();
    let log_range = max_period.ln() - log_min;
    if density_g_cm3.is_none() {
        let duration = minimum_duration(min_period, cadence, None);
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
        let duration = minimum_duration(period, cadence, density_g_cm3);
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
    fn adaptive_grid_preserves_log_spacing_without_density_priors() {
        let grid = adaptive_period_grid(
            0.5,
            20.0,
            1000.0,
            2.0 / 1440.0,
            None,
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
            Some(1.0),
            TEST_BUDGET,
            &mut Vec::new(),
        );
        assert!(grid.len() > 10);
        let first_log_step = (grid[1] / grid[0]).ln();
        let last = grid.len() - 1;
        let last_log_step = (grid[last] / grid[last - 1]).ln();
        assert!(last_log_step > first_log_step);
        for pair in grid.windows(2) {
            let duration = minimum_duration(pair[0], cadence, Some(1.0));
            let timing_drift = span * (pair[1] / pair[0] - 1.0);
            assert!(
                timing_drift <= duration * 0.51,
                "{pair:?}: {timing_drift} > {duration}"
            );
        }
    }
}

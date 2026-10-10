pub(crate) fn median<I: IntoIterator<Item = f64>>(iter: I) -> f64 {
    let mut values: Vec<f64> = iter.into_iter().filter(|v| v.is_finite()).collect();
    if values.is_empty() {
        return 0.0;
    }
    let mid = values.len() / 2;
    let even = values.len().is_multiple_of(2);
    let (lower, middle, _) = values.select_nth_unstable_by(mid, f64::total_cmp);
    if even {
        let left = lower.iter().copied().max_by(f64::total_cmp).unwrap();
        0.5 * left + 0.5 * *middle
    } else {
        *middle
    }
}

pub(crate) fn robust_sigma(values: &[f64]) -> f64 {
    let center = median(values.iter().copied());
    median(values.iter().map(|v| (v - center).abs())) * 1.4826
}

pub(crate) fn phase_distance(a: f64, b: f64, period: f64) -> f64 {
    ((a - b + period / 2.0).rem_euclid(period) - period / 2.0).abs()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn median_handles_even_odd_nonfinite_and_large_values() {
        assert_eq!(median([3.0, 1.0, 2.0]), 2.0);
        assert_eq!(median([4.0, 1.0, 3.0, 2.0]), 2.5);
        assert_eq!(median([f64::NAN, 2.0, f64::INFINITY]), 2.0);
        assert_eq!(median([f64::MAX, f64::MAX]), f64::MAX);
    }
}

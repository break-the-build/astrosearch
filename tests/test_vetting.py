"""Offline tests for the CTOI vetting toolkit (no network, no astrosearch binary)."""
import math
import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from vetting import aliases, diagnostics, fap, gate, stellar  # noqa: E402
from vetting.io import LightCurve, parse_event_epochs  # noqa: E402


def synthetic(period=3.0, epoch=1.0, duration=0.1, depth=0.002, sigma=2e-4, span=54.0,
              cadence=0.01, seed=1, odd_depth=None, secondary=None, parts=None):
    rng = np.random.default_rng(seed)
    time = np.arange(0.0, span, cadence)
    flux = 1.0 + sigma * rng.standard_normal(time.size)
    offset = (time - epoch + period / 2) % period - period / 2
    cycle = np.rint((time - epoch) / period).astype(int)
    inside = np.abs(offset) <= duration / 2
    this_depth = np.where((cycle % 2 == 1) & (odd_depth is not None), odd_depth or 0, depth)
    flux[inside] -= this_depth[inside]
    if secondary:
        phase, sec_depth = secondary
        off2 = (time - epoch - phase * period + period / 2) % period - period / 2
        flux[np.abs(off2) <= duration / 2] -= sec_depth
    part = np.zeros(time.size, dtype=int) if parts is None else parts(time)
    metas = tuple({"sector": i + 1, "author": "SPOC", "crowding_corrected": True,
                   "momentum_dump_times": [], "scattered_light_intervals": []}
                  for i in range(int(part.max()) + 1))
    return LightCurve(time, flux, np.full(time.size, sigma), part, metas)


class DiagnosticsTests(unittest.TestCase):
    def test_event_depths_recover_injected_depth(self):
        lc = synthetic()
        summary = diagnostics.combine(diagnostics.event_depths(lc, 1.0, 3.0, 0.1))
        self.assertAlmostEqual(summary["depth_ppm"], 2000, delta=5 * summary["sigma_ppm"])
        self.assertGreaterEqual(summary["events"], 17)
        self.assertGreater(summary["p_value"], 0.001)

    def test_odd_even_flags_alternating_depths(self):
        self.assertEqual(diagnostics.odd_even(
            diagnostics.event_depths(synthetic(), 1.0, 3.0, 0.1))["status"], "PASS")
        eb = synthetic(odd_depth=0.001)
        self.assertEqual(diagnostics.odd_even(
            diagnostics.event_depths(eb, 1.0, 3.0, 0.1))["status"], "FAIL")

    def test_secondary_scan_finds_eccentric_secondary(self):
        clean = diagnostics.secondary_scan(synthetic(), 1.0, 3.0, 0.1)
        self.assertEqual(clean["status"], "PASS")
        eccentric = diagnostics.secondary_scan(synthetic(secondary=(0.37, 0.0008)), 1.0, 3.0, 0.1)
        self.assertEqual(eccentric["status"], "FAIL")
        self.assertAlmostEqual(eccentric["strongest"]["phase"], 0.37, delta=0.02)

    def test_part_consistency_detects_depth_change_between_sectors(self):
        lc = synthetic(parts=lambda t: (t > 27).astype(int))
        lc.flux[(lc.part == 1) & (np.abs((lc.time - 1 + 1.5) % 3 - 1.5) <= 0.05)] += 0.001
        result = diagnostics.part_consistency(lc, diagnostics.event_depths(lc, 1.0, 3.0, 0.1))
        self.assertEqual(result["status"], "FAIL")
        same = synthetic(parts=lambda t: (t > 27).astype(int))
        self.assertEqual(diagnostics.part_consistency(
            same, diagnostics.event_depths(same, 1.0, 3.0, 0.1))["status"], "PASS")

    def test_quality_overlap_and_grazing(self):
        lc = synthetic()
        events = diagnostics.event_depths(lc, 1.0, 3.0, 0.1)
        flagged = LightCurve(lc.time, lc.flux, lc.error, lc.part,
                             ({**lc.parts[0], "momentum_dump_times": [e["midpoint"] for e in events]},))
        self.assertEqual(diagnostics.quality_overlap(flagged, events, 0.1)["status"], "FAIL")
        self.assertEqual(diagnostics.quality_overlap(lc, events, 0.1)["status"], "PASS")
        self.assertEqual(diagnostics.grazing({"impact": 0.95, "rp_rs": 0.1})["status"], "FAIL")
        self.assertEqual(diagnostics.grazing({"impact": 0.3, "rp_rs": 0.02})["status"], "PASS")


class AliasTests(unittest.TestCase):
    def test_duo_alias_ruled_out_by_flat_data_at_predicted_epoch(self):
        # Transits at 5 and 25; the k=2 alias predicts one at 15, where data are flat.
        lc = synthetic(period=20.0, epoch=5.0, duration=0.2, span=30.0)
        result = {a["k"]: a for a in aliases.duo_aliases(lc, 5.0, 25.0, 0.2, 2000, 60, min_period=3.0)}
        self.assertTrue(result[1]["allowed"])
        self.assertFalse(result[2]["allowed"])

    def test_predicted_windows_propagate_covariance(self):
        windows = aliases.predicted_windows(100.0, 2.0, 115.0, 121.0, 0.001, 0.0001, 1e-8, 0.1)
        self.assertEqual([w["epoch_number"] for w in windows], [8, 9, 10])
        expected = math.sqrt(0.001 ** 2 + 100 * 0.0001 ** 2 + 20 * 1e-8)
        self.assertAlmostEqual(windows[-1]["sigma_days"], expected, places=12)


class StellarTests(unittest.TestCase):
    def test_sun_earth_values(self):
        self.assertAlmostEqual(stellar.transit_density_g_cm3(215.03, 365.256), 1.41, delta=0.02)
        sun = {"rad": 1.0, "e_rad": 0.0, "mass": 1.0, "e_mass": 0.0, "Teff": 5772, "e_Teff": 0,
               "rho": 1.0, "e_rho": 0.0}
        fit = {"rp_rs": 1 / 109.076, "rp_rs_sigma": 0, "a_rs": 215.03, "a_rs_sigma": 0,
               "period": 365.256, "period_sigma": 0}
        planet = stellar.derive_planet(sun, fit)
        self.assertAlmostEqual(planet["radius_earth"]["median"], 1.0, places=3)
        self.assertAlmostEqual(planet["semi_major_axis_au"]["median"], 1.0, delta=0.002)
        self.assertAlmostEqual(planet["insolation_earth"]["median"], 1.0, delta=0.01)
        self.assertAlmostEqual(planet["equilibrium_temperature_k"]["median"], 255, delta=2)
        diluted = stellar.derive_planet(sun, fit, depth_correction=1.21)
        self.assertAlmostEqual(diluted["radius_earth"]["median"], 1.1, places=3)

    def test_density_consistency_statuses(self):
        star = {"rho": 1.0, "e_rho": 0.1}
        consistent = {"median": 1.45, "minus": 0.1, "plus": 0.1}
        self.assertEqual(stellar.density_consistency(consistent, star)["status"], "CONSISTENT")
        tension = {"median": 4.0, "minus": 0.1, "plus": 0.1}
        self.assertEqual(stellar.density_consistency(tension, star)["status"], "TENSION")
        giant = {"median": 0.01, "minus": 0.001, "plus": 0.001}
        self.assertEqual(stellar.density_consistency(giant, star)["status"], "INCONSISTENT")
        self.assertEqual(stellar.density_consistency(None, star)["status"], "NOT_EVALUATED")

    def test_blend_capability_uses_flux_ratio(self):
        star = {"Tmag": 10.0}
        neighbors = [{"ID": 1, "Tmag": 15.0, "dstArcSec": 30.0},   # ~0.01 max depth
                     {"ID": 2, "Tmag": 12.0, "dstArcSec": 100.0}]  # outside 2 pixels
        blends = stellar.blend_capability(star, neighbors, depth_ppm=500)
        self.assertEqual(len(blends), 1)
        self.assertTrue(blends[0]["can_produce_depth"])
        self.assertFalse(stellar.blend_capability(star, neighbors, depth_ppm=20000)[0]["can_produce_depth"])


class FapTests(unittest.TestCase):
    def test_nulls_preserve_noise_and_remove_signal(self):
        lc = synthetic()
        _, inv, _ = fap.inverted(lc)
        self.assertAlmostEqual(np.median(inv), np.median(lc.flux), places=9)
        rng = np.random.default_rng(3)
        time, flux, _ = fap.block_shuffled(lc, rng, 1.0, 3.0, 0.1)
        kept = np.abs((lc.time - 1 + 1.5) % 3 - 1.5) > 0.1
        np.testing.assert_array_equal(time, lc.time[kept])
        np.testing.assert_array_equal(np.sort(flux), np.sort(lc.flux[kept]))
        self.assertFalse(np.array_equal(flux, lc.flux[kept]))

    def test_probability_helpers(self):
        self.assertAlmostEqual(fap.clopper_pearson_upper(0, 300), 1 - 0.05 ** (1 / 300), places=6)
        maxima = np.random.default_rng(1).gumbel(7, 0.5, 400)
        self.assertGreater(fap.gumbel_tail(maxima, 7.0), fap.gumbel_tail(maxima, 9.0))
        self.assertLess(fap.gumbel_tail(maxima, 11.0), 1e-3)


class GateTests(unittest.TestCase):
    def test_verdict_requires_every_automated_criterion(self):
        passing = {"a": {"status": "PASS"}, "b": {"status": "ACTION_REQUIRED"}}
        self.assertTrue(gate.verdict(passing)["ready_for_human_review"])
        blocked = {"a": {"status": "PASS"}, "c": {"status": "NOT_EVALUATED"}}
        decision = gate.verdict(blocked)
        self.assertFalse(decision["ready_for_human_review"])
        self.assertEqual(decision["blockers"], ["c"])

    def test_significance_needs_upper_bound_below_one_percent(self):
        weak = gate.criterion_significance({"fap_upper_95": 0.03}, {"status": "PASS"})
        self.assertEqual(weak["status"], "FAIL")
        strong = gate.criterion_significance({"fap_upper_95": 0.009}, {"status": "PASS"})
        self.assertEqual(strong["status"], "PASS")

    def test_duo_odd_even_is_not_a_blocker(self):
        result = gate.criterion_false_positive(
            {"status": "PASS"}, {"status": "NOT_EVALUATED"}, {"status": "PASS"},
            {"status": "PASS"}, [], 1.0, [{"centroid_separation_arcsec": 5.0}], duo=True)
        self.assertEqual(result["status"], "PASS")
        offset = gate.criterion_false_positive(
            {"status": "PASS"}, {"status": "PASS"}, {"status": "PASS"}, {"status": "PASS"},
            [], 1.0, [{"centroid_separation_arcsec": 40.0}], duo=False)
        self.assertEqual(offset["checks"]["no_centroid_offset"], "FAIL")

    def test_event_evidence_parsing(self):
        self.assertEqual(parse_event_epochs("1468.4:542:4.3:3;1471.6:445:4.1:4"), (1468.4, 1471.6))
        self.assertEqual(parse_event_epochs(""), ())


class RefineTests(unittest.TestCase):
    def test_transit_fit_delta_bic_matches_signal_strength(self):
        from types import SimpleNamespace

        import refine_candidates as refine

        lc = synthetic(depth=0.002, sigma=4e-4, span=24.0, cadence=0.02)
        n = len(lc.time)
        args = SimpleNamespace(limb_darkening=(0.4, 0.2), integration_samples=3, initial_impact=0.4,
                               epoch_window_durations=1.0, fit_window_durations=2.5,
                               period_window_fraction=0.01, robust_scale=1.0, max_evaluations=3000,
                               mcmc_steps=0, mcmc_seed=1)
        row = {"period_days": "3.0", "epoch": "1.0", "duration_days": "0.1", "depth_ppm": "2000",
               "class_rank": "1"}
        fit = refine.fit_candidate(lc.time, lc.flux, lc.error, np.full(n, 0.02),
                                   np.zeros(n, dtype=object).astype(str), row, args)
        self.assertAlmostEqual(fit["fitted_period_days"], 3.0, delta=0.004)
        self.assertAlmostEqual(fit["fitted_radius_ratio"] ** 2, 0.002, delta=0.0008)
        # Expected improvement ~ (depth / sigma)^2 * in-transit points; the old
        # bug mixed transit parameters into the null baselines (~1e14).
        in_transit = 8 * 5
        expected = (0.002 / 4e-4) ** 2 * in_transit
        self.assertLess(fit["delta_bic_null_minus_transit"], 3 * expected)
        self.assertGreater(fit["delta_bic_null_minus_transit"], expected / 3)


if __name__ == "__main__":
    unittest.main()

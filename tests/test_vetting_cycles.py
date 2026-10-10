"""Offline tests for the review-board improvement cycles (no network, no binary).

Each test pins one verified defect found by the known-systems benchmark or
the review, so a regression reintroducing it fails here.
"""
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from vetting import diagnostics, gate, pixel, siblings, stellar, tier  # noqa: E402
from vetting.io import Candidate, LightCurve, normalize_time, read_lightcurve  # noqa: E402


def lightcurve(time, flux, sigma):
    meta = ({"sector": 1, "author": "SPOC", "crowding_corrected": True,
             "momentum_dump_times": [], "scattered_light_intervals": []},)
    return LightCurve(time, flux, np.full(time.size, sigma), np.zeros(time.size, dtype=int), meta)


def box(time, period, epoch, duration, depth):
    offset = (time - epoch + period / 2) % period - period / 2
    return np.where(np.abs(offset) <= duration / 2, depth, 0.0)


def candidate(rank, period, epoch, duration, score, kind="periodic_dip"):
    return Candidate(rank=rank, kind=kind, score=score, period=period, epoch=epoch,
                     duration=duration, depth_ppm=1000.0, n_events=10, event_epochs=(), row={})


class GaiaLookupTests(unittest.TestCase):
    """Cycle 1: a J2000 TIC position misses high-proper-motion stars in Gaia
    DR3 (epoch 2016.0); LHS 3844 moved ~20 arcsec, outside the 2 arcsec cone."""

    def test_query_joins_dr2_neighbourhood_by_source_id(self):
        adql = stellar.gaia_query_by_dr2_id("6385548541499112448")
        self.assertIn("gaiadr3.dr2_neighbourhood", adql)
        self.assertIn("6385548541499112448", adql)
        self.assertIn("ruwe", adql)

    def test_dr2_id_must_be_numeric(self):
        with self.assertRaises(ValueError):
            stellar.gaia_query_by_dr2_id("1 OR 1=1")

    def test_proper_motion_propagation(self):
        # 1000 mas/yr east at dec 60 for 16 years: 16 arcsec on sky, 32 arcsec in RA.
        ra, dec = stellar.propagate_position(10.0, 60.0, 1000.0, 0.0, 2000.0, 2016.0)
        self.assertAlmostEqual((ra - 10.0) * 3600 * math.cos(math.radians(60.0)), 16.0, places=3)
        self.assertAlmostEqual(dec, 60.0, places=9)
        ra, dec = stellar.propagate_position(10.0, 0.0, None, None, 2000.0, 2016.0)
        self.assertEqual((ra, dec), (10.0, 0.0))


class NeighborEpochTests(unittest.TestCase):
    """Cycle 2: neighbor separations use positions at the TESS epoch."""

    def test_separation_moves_with_target_proper_motion(self):
        star = {"ra": 0.0, "dec": 0.0, "pmRA": 0.0, "pmDEC": 2000.0, "Tmag": 10.0}
        neighbor = {"ID": "2", "ra": 0.0, "dec": 40 / 3600, "pmRA": 0.0, "pmDEC": 0.0,
                    "Tmag": 12.0, "dstArcSec": 40.0}
        at_2019 = stellar.neighbors_at_epoch(star, [neighbor], 2019.0)
        self.assertAlmostEqual(at_2019[0]["dstArcSec"], 2.0, delta=0.01)
        self.assertAlmostEqual(at_2019[0]["dstArcSec_j2000"], 40.0, places=6)

    def test_btjd_to_julian_year(self):
        self.assertAlmostEqual(stellar.btjd_to_year(0.0), 2014.9, delta=0.05)
        self.assertAlmostEqual(stellar.btjd_to_year(1325.0), 2018.53, delta=0.05)


class BlendShapeTests(unittest.TestCase):
    """Cycle 3: Vanderburg et al. (2019, ApJL 881, L19) — the ingress/egress
    to total-duration ratio bounds how diluted a transit can be."""

    def test_flat_bottomed_transit_limits_delta_mag(self):
        fit = {"rp_rs": 0.062, "rp_rs_sigma": 0.001, "impact": 0.02, "impact_sigma": 0.05}
        limit = stellar.shape_delta_mag_limit(fit, depth_ppm=3700.0)
        self.assertLess(limit["delta_mag_max"], 1.0)
        self.assertGreater(limit["delta_mag_max"], 0.0)

    def test_grazing_geometry_falls_back_to_full_eclipse_limit(self):
        fit = {"rp_rs": 0.1, "impact": 0.95}
        limit = stellar.shape_delta_mag_limit(fit, depth_ppm=1000.0)
        self.assertAlmostEqual(limit["delta_mag_max"], 2.5 * math.log10(1 / 1e-3), places=6)

    def test_missing_fit_gives_no_limit(self):
        self.assertIsNone(stellar.shape_delta_mag_limit({}, depth_ppm=1000.0)["delta_mag_max"])

    def test_faint_neighbor_excluded_by_shape(self):
        star = {"Tmag": 11.92}
        neighbors = [{"ID": "n", "Tmag": 16.67, "dstArcSec": 41.0}]
        flux_only = stellar.blend_capability(star, neighbors, 3700.0)
        self.assertTrue(flux_only[0]["can_produce_depth"])
        shaped = stellar.blend_capability(star, neighbors, 3700.0, delta_mag_max=0.3)
        self.assertFalse(shaped[0]["can_produce_depth"])
        self.assertTrue(shaped[0]["flux_allows_depth"])
        self.assertFalse(shaped[0]["shape_allows_dilution"])


class SiblingMaskTests(unittest.TestCase):
    """Cycle 4: in multi-planet systems another planet's transits land near
    fixed phases of the candidate and fail its secondary and odd/even tests.
    The geometry mirrors L 98-59 d (7.45 d) and c (3.69 d): a near-2:1 pair."""

    target = candidate(1, 7.45, 1.0, 0.1, 40.0)
    sibling = candidate(2, 3.69, 1.3, 0.1, 30.0)

    def setUp(self):
        rng = np.random.default_rng(4)
        self.time = np.arange(0.0, 60.0, 0.01)
        flux = 1.0 + 2e-4 * rng.standard_normal(self.time.size)
        flux -= box(self.time, 7.45, 1.0, 0.1, 0.002)
        flux -= box(self.time, 3.69, 1.3, 0.1, 0.0015)
        self.lc = lightcurve(self.time, flux, 2e-4)

    def test_unmasked_sibling_fails_odd_even(self):
        # The near-2:1 sibling lands beside alternate transits of the candidate.
        events = diagnostics.event_depths(self.lc, 1.0, 7.45, 0.1)
        self.assertEqual(diagnostics.odd_even(events)["status"], "FAIL")

    def test_masked_sibling_passes_secondary_and_odd_even(self):
        masked, record = siblings.mask_siblings(self.lc, self.target, [self.sibling], baseline=60.0)
        self.assertEqual(len(record["masked"]), 1)
        self.assertEqual(diagnostics.secondary_scan(masked, 1.0, 7.45, 0.1)["status"], "PASS")
        events = diagnostics.event_depths(masked, 1.0, 7.45, 0.1)
        self.assertEqual(diagnostics.odd_even(events)["status"], "PASS")

    def test_same_period_sibling_is_never_masked(self):
        # A same-period dip at phase 0.5 is the eclipsing-binary secondary itself.
        secondary = [candidate(2, 7.4502, 4.725, 0.1, 25.0)]
        _, record = siblings.mask_siblings(self.lc, self.target, secondary, baseline=60.0)
        self.assertEqual(record["masked"], [])
        self.assertEqual(record["kept"][0]["reason"], "commensurate period")

    def test_low_score_siblings_are_not_masked(self):
        weak = [candidate(2, 3.69, 1.3, 0.1, 8.0)]
        _, record = siblings.mask_siblings(self.lc, self.target, weak, baseline=60.0)
        self.assertEqual(record["masked"], [])


class SecondaryCalibrationTests(unittest.TestCase):
    """Cycle 5: the phase scan must account for red noise and the number of
    phases tried (look-elsewhere), as the Kepler Robovetter model-shift test does."""

    def red_noise_curve(self, seed):
        rng = np.random.default_rng(seed)
        time = np.arange(0.0, 54.0, 0.01)
        white = 2e-4 * rng.standard_normal(time.size)
        kernel = np.exp(-np.arange(0, 0.3, 0.01) / 0.08)
        red = np.convolve(rng.standard_normal(time.size), kernel / kernel.sum(), "same") * 6e-4
        flux = 1.0 + white + red - box(time, 3.0, 1.0, 0.1, 0.003)
        return lightcurve(time, flux, 2e-4)

    def test_red_noise_alone_rarely_fails_calibrated_scan(self):
        fails = sum(diagnostics.secondary_scan(self.red_noise_curve(s), 1.0, 3.0, 0.1)["status"] == "FAIL"
                    for s in range(10))
        self.assertLessEqual(fails, 1)

    def test_scan_reports_red_noise_factor_and_threshold(self):
        result = diagnostics.secondary_scan(self.red_noise_curve(0), 1.0, 3.0, 0.1)
        self.assertGreater(result["red_noise_factor"], 1.0)
        self.assertGreater(result["threshold_sigma"], 3.0)

    def test_strong_secondary_still_fails(self):
        lc = self.red_noise_curve(1)
        lc.flux[box(lc.time, 3.0, 2.5, 0.1, 1.0) > 0] -= 0.0015
        self.assertEqual(diagnostics.secondary_scan(lc, 1.0, 3.0, 0.1)["status"], "FAIL")


class PixelCentroidTests(unittest.TestCase):
    """Issue 8: difference-image centroids use significant pixels and an
    event-bootstrap uncertainty; a blended eclipsing neighbor is flagged."""

    def cube(self, eclipse_at, depth=0.01, seed=3):
        rng = np.random.default_rng(seed)
        time = np.arange(0.0, 27.0, 0.02)
        yy, xx = np.indices((11, 11))

        def psf(x0, y0, flux):
            return flux * np.exp(-((xx - x0) ** 2 + (yy - y0) ** 2) / (2 * 0.8 ** 2))

        target, neighbor = psf(5, 5, 1e5), psf(8, 5, 3e4)
        inside = np.abs((time - 1.0 + 1.5) % 3.0 - 1.5) <= 0.05
        frames = np.repeat((target + neighbor + 100.0)[None], time.size, axis=0)
        source = target if eclipse_at == "target" else neighbor
        frames[inside] -= depth * source
        frames += rng.normal(0, 20.0, frames.shape)
        return time, frames

    def test_on_target_transit_centroid_at_target(self):
        time, cube = self.cube("target")
        result = pixel.centroid_with_uncertainty(time, cube, 1.0, 3.0, 0.1)
        self.assertEqual(result["status"], "OK")
        offset = pixel.offset_significance(result["centroid_xy"], result["centroid_sigma_xy"], (5, 5))
        self.assertLess(offset["offset_pixels"], 0.3)

    def test_blended_neighbor_is_significant_offset(self):
        time, cube = self.cube("neighbor", depth=0.05)
        result = pixel.centroid_with_uncertainty(time, cube, 1.0, 3.0, 0.1)
        offset = pixel.offset_significance(result["centroid_xy"], result["centroid_sigma_xy"], (5, 5))
        self.assertGreater(offset["offset_pixels"], 2.0)
        report = {"centroid_separation_arcsec": offset["offset_arcsec"],
                  "centroid_offset_sigma": offset["offset_sigma"]}
        self.assertEqual(pixel.centroid_status(report), "FAIL")

    def test_large_but_insignificant_offset_passes(self):
        self.assertEqual(pixel.centroid_status({"centroid_separation_arcsec": 30.0,
                                                "centroid_offset_sigma": 1.2}), "PASS")
        self.assertEqual(pixel.centroid_status({}), "NOT_EVALUATED")


class TimeSystemTests(unittest.TestCase):
    """Cycle 6: full BJD inputs are normalized to BTJD (BJD - 2457000)."""

    def test_normalize_time(self):
        out, offset = normalize_time(np.array([2458400.5, 2458401.5]))
        self.assertEqual(offset, 2457000.0)
        np.testing.assert_allclose(out, [1400.5, 1401.5])
        out, offset = normalize_time(np.array([1400.5]))
        self.assertEqual(offset, 0.0)

    def test_read_lightcurve_accepts_bjd(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "lc.csv"
            path.write_text("time,flux,flux_err\n2458400.5,1.0,0.001\n2458400.6,1.0,0.001\n")
            lc = read_lightcurve(path)
        np.testing.assert_allclose(lc.time, [1400.5, 1400.6])


class DetectionTierTests(unittest.TestCase):
    """Cycle 7: significance requires the calibrated detection tier, not the
    screening threshold of 6 that 97% of inverted null controls exceed."""

    calibration = {"threshold": 37.83, "target_false_alarm_rate": 0.05,
                   "context": {"transit_search_passes": 1, "flux_origin": "QLP"}}

    def test_tier_check(self):
        self.assertEqual(gate.detection_tier(40.0, self.calibration)["status"], gate.PASS)
        self.assertEqual(gate.detection_tier(20.0, self.calibration)["status"], gate.FAIL)
        self.assertEqual(gate.detection_tier(20.0, None)["status"], gate.NE)

    def test_significance_includes_tier(self):
        fap = {"fap_upper_95": 0.001}
        tier = gate.detection_tier(20.0, self.calibration)
        result = gate.criterion_significance(fap, {"status": gate.PASS}, tier=tier)
        self.assertEqual(result["status"], gate.FAIL)
        self.assertEqual(result["checks"]["score_at_or_above_calibrated_detection_tier"], gate.FAIL)

    def test_jackknife_tier_discounts_single_event_candidates(self):
        calibration = {"statistic": "jackknife_score", "threshold": 23.3}
        even = {"score": "30", "event_evidence": ";".join("1.0:900:5.0:10" for _ in range(9))}
        single = {"score": "30", "event_evidence": "1.0:900:20.0:10;4.0:10:1.0:10;7.0:10:1.0:10"}
        self.assertEqual(gate.detection_tier(even, calibration)["status"], gate.PASS)
        self.assertEqual(gate.detection_tier(single, calibration)["status"], gate.FAIL)
        self.assertEqual(gate.detection_tier(30.0, calibration)["status"], gate.NE)
        self.assertAlmostEqual(tier.jackknife_score(even), 30 * math.sqrt(8 / 9))
        self.assertEqual(tier.jackknife_score({**even, "jackknife_score": "12.5"}), 12.5)

    def test_protocol_never_completes_without_humans(self):
        report = {"provenance": {"git": {"commit": "abc", "dirty": False}},
                  "criteria": {"5_significance": {"detection_tier": {"status": gate.PASS}},
                               "3_no_false_positive_indicators": {"checks": {"no_centroid_offset": gate.PASS}}},
                  "verdict": {"ready_for_human_review": True}, "bls": {"status": gate.PASS}}
        steps = gate.protocol_checklist(report)
        self.assertTrue(all(s["status"] == gate.PASS for s in steps[:5]))
        self.assertTrue(all(s["status"] == gate.ACTION for s in steps[5:]))
        dirty = gate.protocol_checklist({**report, "provenance": {"git": {"commit": "abc", "dirty": True}}})
        self.assertEqual(dirty[0]["status"], gate.FAIL)

    def test_calibration_file_round_trip(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "c.json"
            path.write_text(json.dumps(self.calibration))
            self.assertEqual(gate.load_calibration(path)["threshold"], 37.83)
            self.assertIsNone(gate.load_calibration(Path(tmp) / "missing.json"))


if __name__ == "__main__":
    unittest.main()

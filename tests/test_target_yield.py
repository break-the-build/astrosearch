"""Expected-yield target ranking: physical trends, calibration handling, CLI."""
import csv
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import target_yield as ty  # noqa: E402

M_DWARF = {"tic_id": "1", "radius_solar": "0.4", "teff_k": "3400", "tmag": "12"}


def expected(row, **settings):
    return ty.expected_new_planets(row, MODEL, ty.Settings(**settings))["expected_new_planets"]


MODEL = ty.load_model()


class OccurrenceTests(unittest.TestCase):
    def test_sag13_normalisation_and_breaks(self):
        rate = ty.occurrence(np.array([365.25]), np.array([1.0]), None)
        self.assertAlmostEqual(float(rate[0]), 0.38)
        inner = ty.occurrence(np.array([1.0, 10.0]), np.array([2.0, 2.0]), None)
        self.assertAlmostEqual(float(inner[0] / inner[1]), 10 ** -1.5)
        m_small, g_small = (ty.occurrence(np.array([20.0]), np.array([2.0]), t) for t in (3400, 5700))
        self.assertAlmostEqual(float(m_small[0] / g_small[0]), 2.0)
        m_giant, g_giant = (ty.occurrence(np.array([20.0]), np.array([10.0]), t) for t in (3400, 5700))
        self.assertAlmostEqual(float(m_giant[0] / g_giant[0]), 0.5)


class YieldTests(unittest.TestCase):
    def test_requires_radius_and_noise(self):
        self.assertIsNone(ty.expected_new_planets({"tmag": "10"}, MODEL))
        self.assertIsNone(ty.expected_new_planets({"radius_solar": "0.4"}, MODEL))
        result = ty.expected_new_planets({"radius_solar": "0.4", "cdpp_1h_ppm": "500"}, MODEL)
        self.assertEqual(result["noise_source"], "cdpp_1h_ppm")
        self.assertAlmostEqual(result["probability_new_planet"],
                               1 - math.exp(-result["expected_new_planets"]))

    def test_detectability_falls_with_noise_dilution_and_star_size(self):
        detectable = [ty.expected_new_planets({**M_DWARF, "tmag": str(t)}, MODEL)["expected_detectable_planets"]
                      for t in (9, 11, 13)]
        self.assertGreater(detectable[0], detectable[1])
        self.assertGreater(detectable[1], detectable[2])
        self.assertLess(expected({**M_DWARF, "contamination_ratio": "1"}), expected(M_DWARF))
        self.assertLess(expected({**M_DWARF, "radius_solar": "0.8"}), expected(M_DWARF))

    def test_more_data_helps(self):
        self.assertGreater(expected({**M_DWARF, "sectors": "12"}), expected({**M_DWARF, "sectors": "2"}))
        self.assertEqual(expected({**M_DWARF, "observed_days": "96"}), expected(M_DWARF))

    def test_prior_searches_discount_new_yield(self):
        plain = ty.expected_new_planets(M_DWARF, MODEL)
        spoc = ty.expected_new_planets({**M_DWARF, "spoc_2min": "true"}, MODEL)
        self.assertEqual(plain["expected_detectable_planets"], spoc["expected_detectable_planets"])
        self.assertLess(spoc["expected_new_planets"], plain["expected_new_planets"])
        self.assertLess(plain["expected_new_planets"], plain["expected_detectable_planets"])

    def test_known_transiting_host_is_boosted(self):
        host = ty.expected_new_planets({**M_DWARF, "known_planet_periods": "3.1; 7.2"}, MODEL)
        self.assertTrue(host["known_transiting_host"])
        self.assertGreater(host["expected_new_planets"], 5 * expected(M_DWARF))

    def test_noise_interpolation_extrapolates_monotonically(self):
        sigmas = [ty.noise_at(MODEL, t)[0] for t in (5, 8, 12, 16)]
        self.assertEqual(sigmas, sorted(sigmas))


class CalibrationTests(unittest.TestCase):
    def test_missing_calibration_is_a_clear_error(self):
        with tempfile.TemporaryDirectory() as tmp, self.assertRaises(ValueError) as caught:
            ty.load_model(Path(tmp))
        self.assertIn(ty.CALIBRATION_ENV, str(caught.exception))

    def test_calibration_files_are_consistent(self):
        root = ty.DEFAULT_CALIBRATION
        curve = json.loads((root / "completeness.json").read_text())["curve"]
        fractions = [c["recovered"] / c["n"] for c in curve]
        self.assertEqual(fractions, sorted(fractions))
        self.assertTrue(np.all(np.diff(MODEL.noise_tmag) > 0))


class CliTests(unittest.TestCase):
    def test_ranked_output_and_sky_cells(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            targets = tmp / "targets.csv"
            with targets.open("w", newline="") as stream:
                w = csv.DictWriter(stream, ["tic_id", "radius_solar", "teff_k", "tmag", "ra", "dec"])
                w.writeheader()
                w.writerow({"tic_id": "11", "radius_solar": "0.8", "teff_k": "5000", "tmag": "13", "ra": "10", "dec": "-62"})
                w.writerow({"tic_id": "22", "radius_solar": "0.3", "teff_k": "3300", "tmag": "11", "ra": "12", "dec": "-61"})
                w.writerow({"tic_id": "33", "tmag": "11", "ra": "200", "dec": "40"})
            out, cells = tmp / "ranked.csv", tmp / "cells.csv"
            self.assertEqual(ty.main([str(targets), "--output", str(out), "--sky-cells", str(cells)]), 0)
            ranked = list(csv.DictReader(out.open()))
            self.assertEqual([r["tic_id"] for r in ranked], ["22", "11", "33"])
            self.assertEqual(ranked[2]["expected_new_planets"], "")
            cell_rows = list(csv.DictReader(cells.open()))
            self.assertEqual(len(cell_rows), 1)
            self.assertEqual((cell_rows[0]["targets"], cell_rows[0]["best_tic_id"]), ("2", "22"))

    def test_invalid_arguments(self):
        for extra in (["--min-period", "5", "--max-period", "2"], ["--cell-degrees", "0"]):
            with self.assertRaises(SystemExit):
                ty.main(["missing.csv", "--output", "x.csv", *extra])


if __name__ == "__main__":
    unittest.main()

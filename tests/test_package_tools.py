"""Offline tests for packaged tools that need network or MAST in production:
catalog cross-matching, the pixel-check CLI, transit-window prediction and
the known-systems harness helpers. External services are replaced by fakes."""
import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import benchmark_known_systems as harness  # noqa: E402
import candidate_preflight as preflight  # noqa: E402
import pixel_difference_check  # noqa: E402
from vetting import catalogs  # noqa: E402


class Reply:
    def __init__(self, text):
        self.text = text

    def raise_for_status(self):
        return None


class CatalogTests(unittest.TestCase):
    def test_known_object_check_matches_toi_and_reports_failures(self):
        toi = [{"toi": "100.01", "tid": 4206066, "pl_orbper": 3.18278, "pl_tranmid": 2458468.423,
                "pl_trandurh": 2.0, "tfopwg_disp": "PC"}]
        ctoi_csv = "CTOI,TIC ID,Period (days),Epoch (BJD),Duration (hours)\n4206066.01,4206066,9.9,2458470.0,1.5\n"
        with patch.object(preflight, "tap", side_effect=lambda q: toi if "from toi" in q else []), \
                patch.object(preflight.SESSION, "get", return_value=Reply(ctoi_csv)):
            result = catalogs.known_object_check(4206066, 3.18278, 1468.423, 0.083, [])
        self.assertEqual(result["status"], "FAIL")
        self.assertEqual(len(result["matches"]["toi"]), 1)
        self.assertEqual(result["matches"]["ctoi"], [])

    def test_unknown_candidate_passes_only_when_every_catalog_answered(self):
        with patch.object(preflight, "tap", return_value=[]), \
                patch.object(preflight.SESSION, "get", return_value=Reply("CTOI,TIC ID\n")):
            self.assertEqual(catalogs.known_object_check(1, 3.0, 1400.0, 0.1, [])["status"], "PASS")

        def broken(query):
            raise TimeoutError("TAP stalled")

        with patch.object(preflight, "tap", side_effect=broken), \
                patch.object(preflight.SESSION, "get", return_value=Reply("CTOI,TIC ID\n")):
            result = catalogs.known_object_check(1, 3.0, 1400.0, 0.1, [])
        self.assertEqual(result["status"], "NOT_EVALUATED")
        self.assertIn("toi", result["errors"])

    def test_ctoi_epoch_in_bjd_is_converted(self):
        rows = [{"CTOI": "1.01", "Period (days)": "3.0", "Epoch (BJD)": "2458400.0", "Duration (hours)": "2.4"}]
        matches = catalogs.ctoi_matches(rows, catalogs.candidate_row(3.0, 1400.0, 0.1))
        self.assertEqual(len(matches), 1)


class FakeTpf:
    """A 9x9 stamp with a star at (4, 4) whose transits dim it by 1%."""

    def __init__(self):
        rng = np.random.default_rng(2)
        self.time = SimpleNamespace(value=np.arange(1400.0, 1427.0, 0.02))
        yy, xx = np.indices((9, 9))
        star = 5e4 * np.exp(-((xx - 4) ** 2 + (yy - 4) ** 2) / 1.3)
        cube = np.repeat((star + 50)[None], self.time.value.size, axis=0)
        inside = np.abs((self.time.value - 1401.0 + 1.5) % 3.0 - 1.5) <= 0.05
        cube[inside] -= 0.01 * star
        self.flux = SimpleNamespace(value=cube + rng.normal(0, 10, cube.shape))
        self.quality = np.zeros(self.time.value.size, dtype=int)
        self.wcs = SimpleNamespace(world_to_pixel=lambda coord: (4.0, 4.0),
                                   pixel_to_world=self._to_world)

    @staticmethod
    def _to_world(x, y):
        from astropy.coordinates import SkyCoord
        import astropy.units as u
        return SkyCoord((10 + (x - 4) * 21 / 3600) * u.deg, (20 + (y - 4) * 21 / 3600) * u.deg)


class PixelCliTests(unittest.TestCase):
    def test_run_pixel_check_with_fake_tpf(self):
        search = type("Search", (), {"__len__": lambda self: 1,
                                     "__getitem__": lambda self, i: SimpleNamespace(
                                         download=lambda **k: FakeTpf())})()
        with tempfile.TemporaryDirectory() as tmp, \
                patch("lightkurve.search_targetpixelfile", return_value=search):
            report = pixel_difference_check.run_pixel_check(1, 1, 3.0, 1401.0, 0.1, 10.0, 20.0, Path(tmp),
                                                            pm_ra=0.0, pm_dec=0.0)
            saved = json.loads((Path(tmp) / "pixel_vet.json").read_text())
        self.assertEqual(report["status"], "PASS")
        self.assertLess(saved["centroid_separation_arcsec"], 21.0)

    def test_missing_tpf_is_an_error(self):
        empty = type("Search", (), {"__len__": lambda self: 0})()
        with tempfile.TemporaryDirectory() as tmp, \
                patch("lightkurve.search_targetpixelfile", return_value=empty):
            with self.assertRaises(RuntimeError):
                pixel_difference_check.run_pixel_check(1, 1, 3.0, 1401.0, 0.1, 10.0, 20.0, Path(tmp))


class PredictTransitsTests(unittest.TestCase):
    def run_cli(self, *args):
        return subprocess.run([sys.executable, str(ROOT / "scripts/predict_transits.py"), *args],
                              capture_output=True, text=True)

    def test_windows_widen_with_cycle_count(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "w.csv"
            done = self.run_cli("--period", "3.0", "--epoch", "1400.0", "--period-sigma", "1e-4",
                                "--epoch-sigma", "1e-3", "--after-btjd", "1401", "--through-btjd", "1430",
                                "--output", str(out))
            self.assertEqual(done.returncode, 0, done.stderr)
            rows = list(csv.DictReader(out.read_text().splitlines()))
        self.assertEqual(int(rows[0]["cycle"]), 1)
        self.assertGreater(float(rows[-1]["sigma_days"]), float(rows[0]["sigma_days"]))

    def test_invalid_ephemeris_rejected(self):
        done = self.run_cli("--period", "-1", "--epoch", "1400", "--period-sigma", "0", "--epoch-sigma", "0",
                            "--after-btjd", "1401", "--through-btjd", "1430", "--output", "x.csv")
        self.assertEqual(done.returncode, 2)


class HarnessHelperTests(unittest.TestCase):
    def test_match_rank_and_gate_record(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "c.csv"
            path.write_text("rank,class,period_days,epoch,duration_days,class_rank\n"
                            "1,periodic_variability,5.0,1400.1,0.1,1\n"
                            "2,periodic_dip,4.0,1400.5,0.1,1\n"
                            "3,periodic_dip,3.0,1400.0,0.1,2\n")
            self.assertEqual(harness.match_rank(path, 3.0, 1400.0, 0.1), (2, 3, "periodic_dip"))
            self.assertEqual(harness.match_rank(path, 7.7, 1400.0, 0.1), (None, None, None))
        report = {"criteria": {"3_no_false_positive_indicators": {"status": "PASS",
                                                                  "checks": {"ruwe_below_1_4": "PASS"}},
                               "4_not_already_known": {"status": "FAIL"}},
                  "secondary": {"status": "PASS"}, "odd_even": {"status": "PASS"},
                  "density": {"status": "CONSISTENT"}, "bls": {"status": "PASS"}, "ruwe": 1.1,
                  "siblings": {"masked": [1, 2]}, "products": None}
        record = harness.gate_record(report)
        self.assertEqual(record["fp_ruwe_below_1_4"], "PASS")
        self.assertEqual(record["siblings_masked"], 2)

    def test_run_records_errors_per_target(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest = Path(tmp) / "m.csv"
            manifest.write_text("label,kind,tic,period_days,epoch_btjd,duration_hours,source\n"
                                "X b,planet,1,3.0,1400.0,2.0,test\n")
            with patch.object(harness, "product_specs", return_value=[]):
                harness.run(manifest, Path(tmp) / "out", Path("/nonexistent"), 2, sys.executable, None)
            rows = list(csv.DictReader((Path(tmp) / "out" / "summary.csv").read_text().splitlines()))
        self.assertEqual(rows[0]["status"], "no SPOC 2-min data")


if __name__ == "__main__":
    unittest.main()

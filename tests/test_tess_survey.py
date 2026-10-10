"""Offline end-to-end tests of the TESS batch survey (tess_survey.main).

MAST and Lightkurve are replaced by an in-memory fake that serves synthetic
sector light curves; the real astrosearch binary searches them. The tests
cover target selection, duplicate and ineligible rows, sector selection,
the noisy-sector screen, missing metadata, search timeouts, manifests,
candidate registration, the shared ledger and the cooldown on rerun.
Needs the release binary (ASTROSEARCH_BIN or target/release/astrosearch).
"""
import csv
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import tess_survey  # noqa: E402

BINARY = Path(os.environ.get("ASTROSEARCH_BIN", ROOT / "target/release/astrosearch"))


class FakeTable:
    def __init__(self, rows, columns):
        self.rows, self.colnames = rows, columns

    def __getitem__(self, index):
        return self.rows[index]

    def __len__(self):
        return len(self.rows)


class FakeCurve:
    def __init__(self, time, flux, error):
        self.time = SimpleNamespace(value=time)
        self.flux = SimpleNamespace(value=flux)
        self.flux_err = SimpleNamespace(value=error)

    def remove_nans(self):
        return self

    def normalize(self):
        return self


class FakeSearch:
    def __init__(self, rows, curves, columns):
        self.table = FakeTable(rows, columns)
        self.curves = curves

    def __len__(self):
        return len(self.table)

    def __getitem__(self, selection):
        curve = self.curves[selection[0]]
        return SimpleNamespace(download=lambda **kwargs: curve)


def sector_curve(sector: int, noise: float = 4e-4, seed: int = 0) -> FakeCurve:
    rng = np.random.default_rng(seed + sector)
    start = 1300.0 + 27.4 * sector
    time = start + np.arange(0, 25.0, 600 / 86400)
    flux = 1 + noise * rng.standard_normal(time.size)
    flux[np.abs((time - 1301.1 + 1.6) % 3.2 - 1.6) <= 0.05] -= 0.003
    flux[::97] = np.nan  # raw products carry gaps; the survey must drop them
    return FakeCurve(time, flux, np.full(time.size, noise))


class FakeLightkurve:
    """search_lightcurve keyed by TIC; each product row is one sector."""
    __version__ = "fake"

    def __init__(self, catalog):
        self.catalog = catalog

    def search_lightcurve(self, target, mission=None, author=None):
        tic = target.split()[-1]
        spec = self.catalog.get((tic, author), {"sectors": []})
        columns = spec.get("columns", ["sector", "exptime", "author", "mission"])
        rows, curves = [], []
        for sector, noise in spec["sectors"]:
            rows.append({"sector": sector, "sequence_number": sector, "exptime": 600.0,
                         "author": author, "mission": f"TESS Sector {sector:02d}"})
            curves.append(sector_curve(sector, noise))
        return FakeSearch(rows, curves, columns)


TARGETS = [
    # tic_id, teff, distance, radius, tmag
    ("111", 3400, 20, 0.4, 10.0),   # two good sectors, baseline > minimum: searched
    ("222", 3500, 25, 0.5, 11.0),   # one noisy sector excluded by the screen
    ("333", 3600, 30, 0.5, 11.0),   # only old sectors: completed_no_data
    ("444", 3700, 30, 0.5, 11.0),   # product table lacks sector metadata: retryable error
    ("555", 9000, 30, 0.5, 11.0),   # too hot: skipped
    ("111", 3400, 20, 0.4, 10.0),   # duplicate row
    ("abc", 3400, 20, 0.4, 10.0),   # invalid TIC
]
CATALOG = {
    ("111", "QLP"): {"sectors": [(10, 4e-4), (15, 4e-4)]},
    ("222", "QLP"): {"sectors": [(10, 4e-4), (11, 4e-4), (15, 4e-4), (16, 2e-2)]},
    ("333", "QLP"): {"sectors": [(5, 4e-4), (6, 4e-4)]},
    ("444", "QLP"): {"sectors": [(10, 4e-4), (15, 4e-4)], "columns": ["exptime", "author", "mission"]},
}


@unittest.skipUnless(BINARY.exists(), f"astrosearch binary not found at {BINARY}")
class SurveyEndToEndTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        targets = self.dir / "targets.csv"
        with targets.open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["tic_id", "teff_k", "distance_pc", "radius_solar", "tmag"])
            writer.writerows(TARGETS)
        self.targets = targets
        self.output = self.dir / "out"

    def tearDown(self):
        self.tmp.cleanup()

    def run_survey(self, *extra):
        argv = ["tess_survey.py", "--targets", str(self.targets), "--output", str(self.output),
                "--binary", str(BINARY), "--authors", "QLP", "--min-sector", "12",
                "--min-baseline", "100", "--min-period", "2", "--max-period", "4",
                "--periods", "300", "--threads", "1", "--transit-search-passes", "1",
                "--pause", "0", "--ledger", str(self.dir / "ledger.sqlite"),
                "--cache-dir", str(self.dir / "cache"), *extra]
        fake = FakeLightkurve(CATALOG)
        with patch.object(sys, "argv", argv), \
                patch.object(tess_survey, "dependencies", lambda: (fake, np)), \
                patch.object(tess_survey, "configure_mast_timeout", lambda seconds: None):
            try:
                tess_survey.main()
                code = 0
            except SystemExit as exit_:
                code = exit_.code
        path = self.output / "batch_summary.json"
        if code == 2 or not path.exists():  # argparse rejection: nothing written
            return code, None
        rows = json.loads(path.read_text())
        # Duplicate rows share a TIC; keep the processed entry, not the duplicate marker.
        summary = {str(r["tic_id"]): r for r in rows if r["status"] != "duplicate_tic_id"}
        summary["duplicates"] = [r for r in rows if r["status"] == "duplicate_tic_id"]
        return code, summary

    def test_full_survey_run_and_cooldown(self):
        code, summary = self.run_survey()
        self.assertEqual(code, 1)  # TIC 444 is a retryable error
        self.assertEqual(summary["111"]["status"], "processed")
        self.assertEqual(summary["333"]["status"], "completed_no_data")
        self.assertEqual(summary["444"]["status"], "error")
        self.assertEqual(summary["555"]["status"], "skipped")
        self.assertEqual(summary["abc"]["status"], "error")
        run = summary["111"]["authors"][0]
        self.assertEqual(run["status"], "searched")
        self.assertEqual(run["sectors"], [10, 15])
        self.assertTrue(Path(run["candidates_csv"]).exists())
        with open(run["candidates_csv"], newline="") as stream:
            best = max((r for r in csv.DictReader(stream) if r["class"] == "periodic_dip"),
                       key=lambda r: float(r["score"]))
        self.assertAlmostEqual(float(best["period_days"]), 3.2, delta=0.01)
        manifest = json.loads((self.output / "TIC111" / "manifest.json").read_text())
        self.assertEqual(manifest["search_config"]["min_sector"], 12)
        self.assertTrue((self.output / "TIC111" / "candidate_registry_ingest.json").exists())
        noisy = summary["222"]["authors"][0]
        self.assertEqual(noisy["excluded_noisy_sectors"], [16])
        self.assertNotIn(16, noisy["sectors"])
        with (self.output / "local_decisions.csv").open() as stream:
            decisions = list(csv.DictReader(stream))
        states = {r["input_row"]: r["status"] for r in decisions}
        self.assertEqual(states["6"], "duplicate_tic_id")
        # Default ranking: expected new planets, recorded for every eligible row.
        for row in decisions:
            self.assertEqual(row["expected_new_planets"] != "", row["eligible"] == "True", row)
        self.assertEqual(len(summary["duplicates"]), 1)
        metadata = json.loads((self.output / "run_metadata.json").read_text())
        self.assertEqual(metadata["status"], "partial")
        self.assertEqual(metadata["ranking"], "expected_new_planets")
        self.assertTrue((self.output / "target_ranking.csv").exists())

        # Rerun: completed targets are deferred by the ledger cooldown; the
        # retryable error is attempted again.
        code, summary = self.run_survey()
        self.assertEqual(summary["111"]["status"], "ledger_deferred")
        self.assertEqual(summary["333"]["status"], "ledger_deferred")
        self.assertEqual(summary["444"]["status"], "error")

    def test_search_timeout_is_retryable(self):
        code, summary = self.run_survey("--search-timeout", "0.01", "--max-targets", "1",
                                        "--ranking", "snr")  # selects TIC 111
        self.assertEqual(code, 1)
        run = summary["111"]["authors"][0]
        self.assertEqual((run["status"], run["reason"]), ("search_failed", "search timeout"))
        self.assertEqual(summary["111"]["status"], "error")

    def test_fast_search_profile_finds_the_signal(self):
        code, summary = self.run_survey("--search-profile", "fast", "--max-targets", "1",
                                        "--ranking", "snr", "--periods", "0")
        run = summary["111"]["authors"][0]
        self.assertEqual(run["status"], "searched")
        log = (self.output / "TIC111" / "QLP" / "search.log")
        with open(run["candidates_csv"], newline="") as stream:
            best = max((r for r in csv.DictReader(stream) if r["class"] == "periodic_dip"),
                       key=lambda r: float(r["score"]))
        self.assertAlmostEqual(float(best["period_days"]), 3.2, delta=0.01)
        self.assertIn("two-stage search (fast folding)", log.read_text())

    def test_snr_ranking_remains_available(self):
        self.run_survey("--ranking", "snr", "--search-timeout", "0.01", "--max-targets", "1")
        metadata = json.loads((self.output / "run_metadata.json").read_text())
        self.assertEqual(metadata["ranking"], "local_snr_proxy")
        with (self.output / "local_decisions.csv").open() as stream:
            self.assertTrue(all(r["expected_new_planets"] == "" for r in csv.DictReader(stream)))

    def test_invalid_arguments_are_rejected(self):
        for extra in (["--max-targets", "0"], ["--search-timeout", "-1"],
                      ["--max-sector-noise-ratio", "0.5"], ["--assumed-sectors", "0"],
                      ["--max-period", "1"], ["--calibration-dir", str(self.dir / "missing")],
                      ["--search-profile", "fast"]):  # run_survey passes --periods 300
            code, summary = self.run_survey(*extra)
            self.assertEqual((code, summary), (2, None), extra)


class SurveyHelperTests(unittest.TestCase):
    def test_tic_parsing_and_priority(self):
        self.assertEqual(tess_survey.tic({"tic_id": "TIC 123.0"}), "123")
        with self.assertRaises(ValueError):
            tess_survey.tic({"tic_id": "-5"})
        self.assertGreater(tess_survey.priority({"distance_pc": "10", "radius_solar": "0.2"}, 10, 2000),
                           tess_survey.priority({}, 1, 100))

    def test_snr_proxy_prefers_cdpp_and_applies_contamination(self):
        base = {"radius_solar": "0.5", "tmag": "10"}
        self.assertIsNotNone(tess_survey.local_snr_proxy(base))
        with_cdpp = tess_survey.local_snr_proxy({**base, "cdpp_1h_ppm": "100"})
        diluted = tess_survey.local_snr_proxy({**base, "cdpp_1h_ppm": "100", "contamination_ratio": "1"})
        self.assertAlmostEqual(diluted, with_cdpp / 2)
        self.assertIsNone(tess_survey.local_snr_proxy({"tmag": "10"}))

    def test_sector_noise_density(self):
        rng = np.random.default_rng(1)
        times = np.arange(0, 1, 0.001)
        quiet = tess_survey.sector_noise_density(times, 1 + 1e-4 * rng.standard_normal(times.size), 120)
        loud = tess_survey.sector_noise_density(times, 1 + 1e-3 * rng.standard_normal(times.size), 120)
        self.assertGreater(loud, 5 * quiet)
        self.assertIsNone(tess_survey.sector_noise_density(times[:10], times[:10], 120))

    def test_atomic_writers(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "x.csv"
            tess_survey.atomic_write_csv(path, ["a"], [{"a": 1}])
            self.assertEqual(path.read_text().split(), ["a", "1"])
            tess_survey.atomic_savetxt(Path(tmp) / "v.csv", np.ones((2, 2)), np, delimiter=",")
            self.assertEqual(len(list(Path(tmp).iterdir())), 2)


if __name__ == "__main__":
    unittest.main()

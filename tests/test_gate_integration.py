"""End-to-end: astrosearch search -> ctoi_gate (offline) on a synthetic transit.

Checks the pipeline contract rather than numbers that depend on the noise
draw: the injected signal is found, every criterion and protocol step is
reported with a legal status, provenance and input accounting are present,
nothing claims readiness without the human steps, and gate.json (written
last, atomically) is complete. Needs the release binary; skipped otherwise.
"""
import csv
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BINARY = Path(os.environ.get("ASTROSEARCH_BIN", ROOT / "target/release/astrosearch"))
SEARCH = ["--min-period", "2", "--max-period", "4", "--periods", "400", "--threads", "2",
          "--transit-search-passes", "1"]
STATUSES = {"PASS", "FAIL", "NOT_EVALUATED", "ACTION_REQUIRED"}


@unittest.skipUnless(BINARY.exists(), f"astrosearch binary not found at {BINARY}")
class GateIntegrationTests(unittest.TestCase):
    def test_search_then_offline_gate(self):
        rng = np.random.default_rng(7)
        time = np.arange(1500.0, 1527.0, 2 / 1440)
        flux = 1 + 3e-4 * rng.standard_normal(time.size)
        flux[np.abs((time - 1501.3 + 1.6) % 3.2 - 1.6) <= 0.05] -= 0.003
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            curve = tmp / "lc.csv"
            np.savetxt(curve, np.column_stack([time, flux, np.full(time.size, 3e-4)]),
                       delimiter=",", header="time,flux,flux_err", comments="", fmt="%.8f")
            candidates = tmp / "candidates.csv"
            subprocess.run([str(BINARY), "--input", str(curve), "--output", str(candidates), *SEARCH],
                           check=True, capture_output=True)
            with candidates.open() as stream:
                rows = [r for r in csv.DictReader(stream) if r["class"] == "periodic_dip"]
            best = max(rows, key=lambda r: float(r["score"]))
            self.assertAlmostEqual(float(best["period_days"]), 3.2, delta=0.01)
            self.assertIn(best["tier"], ("detection", "screening"))
            out = tmp / "gate"
            done = subprocess.run(
                [sys.executable, str(ROOT / "scripts/ctoi_gate.py"), "--tic", "1", "--offline",
                 "--lightcurve", str(curve), "--candidates", str(candidates), "--rank", best["rank"],
                 "--binary", str(BINARY), "--search-args", " ".join(SEARCH), "--null-trials", "3",
                 "--workers", "3", "--threads", "1", "--mcmc-steps", "0", "--pixel-sectors", "0",
                 "--out-dir", str(out)], capture_output=True, text=True)
            self.assertEqual(done.returncode, 0, done.stderr[-2000:])
            report = json.loads((out / "gate.json").read_text())
            for name in ("gate.md", "ctoi_row.csv", "ephemeris_windows.csv", "fit.csv"):
                self.assertTrue((out / name).exists(), name)
        self.assertEqual(set(report["criteria"]), {
            "1_multisector_multiproduct", "2_physical_fit_density", "3_no_false_positive_indicators",
            "4_not_already_known", "5_significance", "6_preregistered_ephemeris"})
        for result in report["criteria"].values():
            self.assertIn(result["status"], STATUSES)
        self.assertFalse(report["verdict"]["ready_for_human_review"])  # offline: catalogs unknown
        steps = report["verification_protocol"]
        self.assertEqual([s["step"] for s in steps], list(range(1, 9)))
        self.assertTrue(all(s["status"] == "ACTION_REQUIRED" for s in steps[5:]))
        self.assertEqual(report["fap"]["nulls"], 3)
        self.assertEqual(report["input_rows"]["rows_skipped"], 0)
        self.assertIn("sha256", report["provenance"]["inputs"]["lightcurve"])
        self.assertAlmostEqual(report["fit"]["period"], 3.2, delta=0.01)


if __name__ == "__main__":
    unittest.main()

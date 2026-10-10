"""The known-answer benchmark runner records per-case resource metrics, keeps
the search's provenance next to the candidates (not under the staging name),
and resumes from a matching sidecar. Needs the release binary."""
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


@unittest.skipUnless(BINARY.exists(), f"astrosearch binary not found at {BINARY}")
class BenchmarkRunnerTests(unittest.TestCase):
    def test_metrics_provenance_and_resume(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            rng = np.random.default_rng(3)
            time = np.arange(0.0, 20.0, 0.01)
            np.savetxt(tmp / "lc.csv", np.column_stack([time, 1 + 1e-3 * rng.standard_normal(time.size)]),
                       delimiter=",", header="time,flux", comments="", fmt="%.6f")
            manifest = tmp / "run-manifest.csv"
            with manifest.open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=["case_id", "data_csv", "candidates_csv"])
                writer.writeheader()
                writer.writerow({"case_id": "c1", "data_csv": "lc.csv",
                                 "candidates_csv": "results/c1/candidates.csv"})
            command = [sys.executable, str(ROOT / "scripts/run_known_answer_benchmark.py"),
                       "--manifest", str(manifest), "--binary", str(BINARY), "--threads", "1",
                       "--min-period", "0.5", "--max-period", "2", "--periods", "50",
                       "--transit-search-passes", "1"]
            first = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(first.returncode, 0, first.stderr)
            self.assertIn("searched", first.stdout)
            case = tmp / "results/c1"
            sidecar = json.loads((case / "candidates.run.json").read_text())
            metrics = sidecar["metrics"]
            self.assertGreater(metrics["wall_seconds"], 0)
            self.assertEqual((metrics["exit_code"], metrics["timed_out"]), (0, False))
            self.assertEqual(metrics["observations"], 2000)
            self.assertTrue((case / "candidates.csv.provenance.json").exists())
            self.assertEqual([p.name for p in case.iterdir() if p.name.startswith(".")], [])
            second = subprocess.run(command, capture_output=True, text=True)
            self.assertIn("cached", second.stdout)


if __name__ == "__main__":
    unittest.main()

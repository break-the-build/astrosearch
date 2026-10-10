"""Differential tests of hand-written numerical routines against trusted
references (SciPy, hashlib). References can be wrong too; these pin
agreement on well-defined quantities where both are standard formulas."""
import hashlib
import json
import os
import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import evaluate_benchmark as ev  # noqa: E402
from vetting import fap  # noqa: E402

BINARY = Path(os.environ.get("ASTROSEARCH_BIN", ROOT / "target/release/astrosearch"))


class IntervalTests(unittest.TestCase):
    def test_wilson_matches_scipy(self):
        rng = random.Random(1)
        for _ in range(200):
            n = rng.randint(1, 500)
            k = rng.randint(0, n)
            ours = ev.wilson(k, n)
            ref = stats.binomtest(k, n).proportion_ci(confidence_level=0.95, method="wilson")
            self.assertAlmostEqual(ours[0], ref.low, places=9, msg=(k, n))
            self.assertAlmostEqual(ours[1], ref.high, places=9, msg=(k, n))

    def test_clopper_pearson_upper_matches_scipy_one_sided(self):
        rng = random.Random(2)
        for _ in range(200):
            n = rng.randint(1, 400)
            k = rng.randint(0, n - 1)
            # One-sided 95% upper bound equals the upper end of a two-sided 90% exact interval.
            ref = stats.binomtest(k, n).proportion_ci(confidence_level=0.90, method="exact").high
            self.assertAlmostEqual(fap.clopper_pearson_upper(k, n), ref, places=9, msg=(k, n))
        self.assertEqual(fap.clopper_pearson_upper(5, 5), 1.0)


@unittest.skipUnless(BINARY.exists(), f"astrosearch binary not found at {BINARY}")
class ProvenanceDigestTests(unittest.TestCase):
    def test_rust_sha256_matches_hashlib_on_random_inputs(self):
        rng = random.Random(3)
        with tempfile.TemporaryDirectory() as tmp:
            for case in range(5):
                lines = ["time,flux"] + [f"{i * 0.02 + rng.random() * 1e-3:.6f},{1 + rng.gauss(0, 1e-3):.6f}"
                                         for i in range(rng.randint(40, 4000))]
                # Vary the length mod 64 to cover SHA-256 padding boundaries.
                text = "\n".join(lines) + "\n" + "#" * rng.randint(0, 130) + "\n"
                path = Path(tmp) / f"c{case}.csv"
                path.write_text(text)
                out = Path(tmp) / f"o{case}.csv"
                subprocess.run([str(BINARY), "--input", str(path), "--output", str(out), "--periods", "10",
                                "--durations", "1", "--transit-search-passes", "1", "--threads", "1",
                                "--min-period", "0.2", "--max-period", "0.3"], check=True, capture_output=True)
                record = json.loads(Path(str(out) + ".provenance.json").read_text())
                self.assertEqual(record["input_sha256"], hashlib.sha256(path.read_bytes()).hexdigest())


if __name__ == "__main__":
    unittest.main()

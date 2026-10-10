"""Time budgets for the gate's slow stages: null searches (FAP) and
TRICERATOPS (FPP) must end in a recorded outcome, never a hang or a crash."""
import math
import stat
import sys
import tempfile
import time
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from vetting import fap, fpp, gate  # noqa: E402
from vetting.budget import run_with_deadline  # noqa: E402
from vetting.io import LightCurve  # noqa: E402

FAKE_SEARCH = """#!{python}
import sys, time
time.sleep({sleep})
out = sys.argv[sys.argv.index("--output") + 1]
open(out, "w").write("rank,class,score\\n1,periodic_dip,7.5\\n")
"""


def fake_binary(directory: Path, sleep: float) -> Path:
    path = directory / f"fake_search_{str(sleep).replace('.', '_')}"
    path.write_text(FAKE_SEARCH.format(python=sys.executable, sleep=sleep))
    path.chmod(path.stat().st_mode | stat.S_IEXEC)
    return path


def light_curve() -> LightCurve:
    rng = np.random.default_rng(0)
    time_ = np.arange(0.0, 27.0, 0.02)
    return LightCurve(time_, 1 + 1e-3 * rng.standard_normal(time_.size), np.full(time_.size, 1e-3),
                      np.zeros(time_.size, dtype=int), ())


class DeadlineTests(unittest.TestCase):
    def test_result_timeout_and_error(self):
        self.assertEqual(run_with_deadline(math.sqrt, (4.0,), timeout=60), 2.0)
        start = time.monotonic()
        with self.assertRaises(TimeoutError):
            run_with_deadline(time.sleep, (30,), timeout=1.0)
        self.assertLess(time.monotonic() - start, 15)
        with self.assertRaises(RuntimeError) as caught:
            run_with_deadline(math.sqrt, (-1.0,), timeout=60)
        self.assertIn("ValueError", str(caught.exception))

    def test_fpp_budget_returns_not_evaluated(self):
        result = fpp.triceratops_fpp_budgeted(1, light_curve(), [1], 3.0, 1.0, 0.1, 1000.0,
                                              timeout=0.05)
        self.assertEqual(result["status"], "NOT_EVALUATED")


class FapBudgetTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        self.lc = light_curve()

    def tearDown(self):
        self.tmp.cleanup()

    def run_fap(self, sleep, **kwargs):
        return fap.empirical_fap(self.lc, fake_binary(self.dir, sleep), 10.0, 1.0, 3.0, 0.1, [],
                                 trials=kwargs.pop("trials", 24), workers=8, threads=1,
                                 work_dir=self.dir / f"nulls{sleep}", **kwargs)

    def test_completed_nulls_give_a_fap(self):
        result = self.run_fap(0.0, min_nulls=20)
        self.assertEqual(result["status"], "OK")
        self.assertEqual((result["nulls"], result["exceeding"]), (24, 0))
        self.assertEqual(result["budget"]["failures"], {})

    def test_hanging_null_searches_time_out(self):
        start = time.monotonic()
        result = self.run_fap(30.0, null_timeout=0.5, min_nulls=20)
        self.assertLess(time.monotonic() - start, 20)
        self.assertEqual(result["status"], "NOT_EVALUATED")
        self.assertEqual(result["budget"]["failures"], {"timed out": 24})
        criterion = gate.criterion_significance(result, {"status": gate.PASS})
        self.assertEqual(criterion["checks"]["empirical_fap_upper95_below_1pct"], gate.NE)

    def test_sequential_test_stops_once_clearly_not_significant(self):
        # Every null scores 7.5; a candidate at 5.0 is exceeded by all of them.
        weak = fap.empirical_fap(self.lc, fake_binary(self.dir, 0.0), 5.0, 1.0, 3.0, 0.1, [],
                                 trials=300, workers=8, threads=1, work_dir=self.dir / "weak")
        self.assertEqual(weak["status"], "OK")
        self.assertLess(weak["nulls"], 60)
        self.assertIn("stopped after", weak["budget"]["stopped_early"])
        self.assertGreater(weak["fap_lower_99"], 0.01)
        criterion = gate.criterion_significance(weak, {"status": gate.PASS})
        self.assertEqual(criterion["checks"]["empirical_fap_upper95_below_1pct"], gate.FAIL)
        full = fap.empirical_fap(self.lc, fake_binary(self.dir, 0.0), 5.0, 1.0, 3.0, 0.1, [],
                                 trials=40, workers=8, threads=1, work_dir=self.dir / "full",
                                 sequential=False)
        self.assertEqual(full["nulls"], 40)

    def test_strong_candidate_is_never_stopped_early(self):
        strong = fap.empirical_fap(self.lc, fake_binary(self.dir, 0.0), 50.0, 1.0, 3.0, 0.1, [],
                                   trials=48, workers=8, threads=1, work_dir=self.dir / "strong")
        self.assertEqual((strong["nulls"], strong["exceeding"]), (48, 0))
        self.assertIsNone(strong["budget"]["stopped_early"])

    def test_stage_budget_skips_remaining_nulls(self):
        result = self.run_fap(0.6, budget_seconds=1.5, trials=40, min_nulls=20)
        failures = result["budget"]["failures"]
        self.assertGreater(failures.get("budget exhausted before start", 0), 0)
        self.assertLess(result["budget"]["completed"], 40)
        self.assertEqual(result["status"], "NOT_EVALUATED")


if __name__ == "__main__":
    unittest.main()

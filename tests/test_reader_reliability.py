"""Reader reliability: invalid input must be dropped or rejected, never
averaged into a vetting statistic."""
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from vetting import diagnostics  # noqa: E402
from vetting.io import read_lightcurve, read_table_with_stats  # noqa: E402


def eclipsing_binary(nan_every: int | None = None) -> str:
    rng = np.random.default_rng(1)
    time = np.arange(0.0, 54.0, 0.01)
    flux = 1.0 + 2e-4 * rng.standard_normal(time.size)
    flux[np.abs((time - 1.0 + 1.5) % 3.0 - 1.5) <= 0.05] -= 0.002
    flux[np.abs((time - 2.5 + 1.5) % 3.0 - 1.5) <= 0.05] -= 0.001  # secondary at phase 0.5
    text = ["time,flux,flux_err"]
    for i, (t, f) in enumerate(zip(time, flux)):
        value = "nan" if nan_every and i % nan_every == 0 else f"{f:.8f}"
        text.append(f"{t:.6f},{value},0.0002")
    return "\n".join(text) + "\n"


class ReaderReliabilityTests(unittest.TestCase):
    def read(self, text: str, name="lc.csv"):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / name
            path.write_text(text)
            return read_lightcurve(path)

    def test_nan_rows_do_not_flip_secondary_verdict(self):
        clean = self.read(eclipsing_binary())
        dirty = self.read(eclipsing_binary(nan_every=50))
        self.assertEqual(dirty.stats.rows_skipped, 108)
        self.assertTrue(np.all(np.isfinite(dirty.flux)))
        for lc in (clean, dirty):
            self.assertEqual(diagnostics.secondary_scan(lc, 1.0, 3.0, 0.1)["status"], "FAIL")
            depth = diagnostics.combine(diagnostics.event_depths(lc, 1.0, 3.0, 0.1))["depth_ppm"]
            self.assertTrue(np.isfinite(depth))

    def test_quality_flags_and_spoc_column_names(self):
        lines = ["TIME,PDCSAP_FLUX,PDCSAP_FLUX_ERR,SAP_FLUX,QUALITY"]
        lines += [f"{1000 + i * 0.02:.4f},{1.0:.4f},0.001,5.0,{512 if i % 10 == 0 else 0}"
                  for i in range(100)]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "spoc.csv"
            path.write_text("\n".join(lines) + "\n")
            data, stats = read_table_with_stats(path)
        self.assertEqual(stats.rows_skipped, 10)
        self.assertEqual(len(data), 90)
        np.testing.assert_allclose(data[:, 1], 1.0)  # PDCSAP preferred over SAP
        np.testing.assert_allclose(data[:, 2], 0.001)  # its paired error

    def test_too_few_valid_rows_is_an_error_not_a_curve(self):
        text = "time,flux\n" + "".join(f"{i},nan\n" for i in range(100)) + "1,1\n"
        with self.assertRaises(ValueError):
            self.read(text)

    def test_missing_columns_are_reported(self):
        with self.assertRaises(ValueError):
            self.read("t_mid,counts\n1,2\n")

    def test_duplicate_times_merge_to_median(self):
        lines = ["time,flux,flux_err"] + [f"{i * 0.02:.4f},1.0,0.001" for i in range(40)]
        lines += ["0.0000,3.0,0.001", "0.0000,5.0,0.001"]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "dup.csv"
            path.write_text("\n".join(lines) + "\n")
            data, stats = read_table_with_stats(path)
        self.assertEqual(stats.duplicates_merged, 2)
        self.assertEqual(len(data), 40)
        self.assertEqual(data[0, 1], 3.0)  # median of 1, 3, 5


if __name__ == "__main__":
    unittest.main()

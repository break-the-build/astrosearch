"""Differential fuzz test: the Rust search and the Python vetting reader must
accept, skip and merge exactly the same rows of a light-curve file.

The gate vets the same file the search scanned; if the two readers disagree,
a vetting verdict can rest on samples the search never saw. A NaN-flux row
once flipped the secondary-eclipse verdict of an eclipsing binary to PASS
because only the Rust reader dropped it.

Needs the release binary (ASTROSEARCH_BIN or target/release/astrosearch);
skipped when it is absent.
"""
import json
import os
import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from vetting.io import read_table_with_stats  # noqa: E402

BINARY = Path(os.environ.get("ASTROSEARCH_BIN", ROOT / "target/release/astrosearch"))
FAST = ["--periods", "10", "--durations", "1", "--transit-search-passes", "1", "--threads", "1",
        "--min-period", "0.5", "--max-period", "1.0"]
HEADERS = [
    ("time", "flux", "flux_err", None),
    ("TIME", "PDCSAP_FLUX", "PDCSAP_FLUX_ERR", "QUALITY"),
    ("btjd", "sap_flux", "sap_flux_err", "quality_flag"),
    ("bjd", "relative_flux", "uncertainty", None),
]
BAD_FLUX = ["nan", "NaN", "inf", "-inf", "-1.0", "0", "", "abc", "1e400"]


def fuzz_file(rng: random.Random, path: Path) -> None:
    time_name, flux_name, error_name, quality_name = rng.choice(HEADERS)
    names = [time_name, flux_name, error_name] + ([quality_name] if quality_name else []) + ["label"]
    rng.shuffle(names)
    delimiter = rng.choice([",", "\t"])
    lines = []
    if rng.random() < 0.5:
        lines.append("# exported light curve")
    lines.append(delimiter.join(names))
    t = 1000.0 + rng.random() * 100
    for _ in range(rng.randint(25, 160)):
        t += rng.choice([0.0, 0.01, 0.02, 0.02, 0.02])  # 0.0 makes exact duplicates
        cells = {time_name: f"{t:.6f}", flux_name: f"{1 + rng.gauss(0, 1e-3):.8f}",
                 error_name: f"{abs(rng.gauss(1e-3, 2e-4)):.6f}", "label": rng.choice(["a", "b c"])}
        if quality_name:
            cells[quality_name] = rng.choice(["0", "0", "0", "0.0", "1", "512", "x", ""])
        roll = rng.random()
        if roll < 0.10:
            cells[flux_name] = rng.choice(BAD_FLUX)
        elif roll < 0.15:
            cells[time_name] = rng.choice(["nan", "", "t", "inf"])
        elif roll < 0.20:
            cells[error_name] = rng.choice(["nan", "-1", "0", "inf", ""])
        row = [cells[n] for n in names]
        if rng.random() < 0.05:
            row = row[: rng.randint(1, len(row) - 1)]  # truncated record
        if rng.random() < 0.05:
            row = [f'"{v}"' if "\t" not in v and '"' not in v else v for v in row]
        lines.append(delimiter.join(f" {v}" if rng.random() < 0.05 else v for v in row))
        if rng.random() < 0.03:
            lines.append(rng.choice(["", "# note", "   "]))
    path.write_text(rng.choice(["\n", "\r\n"]).join(lines) + "\n")


def rust_accounting(path: Path):
    out = path.with_suffix(".candidates.csv")
    done = subprocess.run([str(BINARY), "--input", str(path), "--output", str(out), *FAST],
                          capture_output=True, text=True)
    if done.returncode != 0:
        return None
    record = json.loads(Path(str(out) + ".provenance.json").read_text())
    skipped = removed = 0
    for warning in record["warnings"]:
        words = warning.split()
        if warning.startswith("skipped"):
            skipped = int(words[1])
        elif warning.startswith("removed"):
            removed = int(words[1])
    return record["observations"], skipped, removed


@unittest.skipUnless(BINARY.exists(), f"astrosearch binary not found at {BINARY}")
class ParserParityTests(unittest.TestCase):
    def test_rust_and_python_readers_agree_on_fuzzed_files(self):
        rng = random.Random(20261010)
        with tempfile.TemporaryDirectory() as tmp:
            compared = 0
            for case in range(60):
                path = Path(tmp) / f"case{case}.csv"
                fuzz_file(rng, path)
                rust = rust_accounting(path)
                try:
                    data, stats = read_table_with_stats(path)
                except ValueError:
                    data, stats = None, None
                if rust is None:
                    # Rust refused the file (e.g. < 30 valid rows); Python must not
                    # silently produce a usable light curve from it.
                    self.assertTrue(data is None or len(data) < 30, path.read_text()[:400])
                    continue
                observations, skipped, removed = rust
                with self.subTest(case=case):
                    self.assertIsNotNone(data)
                    self.assertEqual(len(data), observations)
                    self.assertEqual(stats.rows_skipped, skipped)
                    self.assertEqual(stats.duplicates_merged, removed)
                compared += 1
            self.assertGreater(compared, 35)


if __name__ == "__main__":
    unittest.main()

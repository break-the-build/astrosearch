"""Mutation fuzzing of the astrosearch CLI input boundary.

Corrupted light-curve files (byte flips, truncation, inserted delimiters and
quotes, invalid UTF-8, binary noise, extreme numbers) must either produce a
normal result or a clean error (exit status 2 with an `astrosearch:` message, src/main.rs).
A signal (release builds abort on panic) or any other exit status is a defect.
Needs the release binary; skipped otherwise.
"""
import os
import random
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BINARY = Path(os.environ.get("ASTROSEARCH_BIN", ROOT / "target/release/astrosearch"))
FAST = ["--periods", "10", "--durations", "1", "--transit-search-passes", "1", "--threads", "1",
        "--min-period", "0.5", "--max-period", "1.0"]


def seed_file(rng: random.Random) -> bytes:
    lines = ["time,flux,flux_err,quality"]
    t = 1400.0
    for _ in range(200):
        t += 0.02
        lines.append(f"{t:.6f},{1 + rng.gauss(0, 1e-3):.8f},0.001,0")
    return ("\n".join(lines) + "\n").encode()


def mutate(rng: random.Random, data: bytes) -> bytes:
    data = bytearray(data)
    for _ in range(rng.randint(1, 6)):
        kind = rng.randrange(8)
        position = rng.randrange(len(data)) if data else 0
        if kind == 0 and data:
            data[position] = rng.randrange(256)
        elif kind == 1:
            data = data[:position]
        elif kind == 2:
            data[position:position] = rng.choice([b",", b'"', b"\t", b"\n", b"\r", b"#"])
        elif kind == 3:
            data[position:position] = bytes(rng.randrange(256) for _ in range(rng.randint(1, 40)))
        elif kind == 4:
            data[position:position] = rng.choice([b"1e308", b"-1e308", b"1e-320", b"nan", b"inf",
                                                  b"99999999999999999999999"])
        elif kind == 5:
            data[position:position] = b"\xff\xfe\xc0\x80"  # invalid UTF-8
        elif kind == 6 and data:
            line = data.find(b"\n", position)
            data[position:position] = data[position:line if line > 0 else len(data)] * 3
        else:
            data = bytearray(b"\xef\xbb\xbf") + data  # BOM
    return bytes(data)


@unittest.skipUnless(BINARY.exists(), f"astrosearch binary not found at {BINARY}")
class CliFuzzTests(unittest.TestCase):
    def test_corrupted_inputs_fail_cleanly(self):
        rng = random.Random(4206066)
        outcomes = {"ok": 0, "error": 0}
        with tempfile.TemporaryDirectory() as tmp:
            for case in range(120):
                path = Path(tmp) / f"f{case}.csv"
                path.write_bytes(mutate(rng, seed_file(rng)) if case else b"")
                done = subprocess.run([str(BINARY), "--input", str(path), "--output",
                                       str(Path(tmp) / f"o{case}.csv"), *FAST],
                                      capture_output=True, timeout=120)
                with self.subTest(case=case):
                    self.assertIn(done.returncode, (0, 2), f"status {done.returncode}: "
                                  f"{done.stderr[-500:]!r}\ninput head {path.read_bytes()[:200]!r}")
                    if done.returncode == 2:
                        self.assertIn(b"astrosearch:", done.stderr)
                outcomes["ok" if done.returncode == 0 else "error"] += 1
        # The fuzzer must exercise both paths to be meaningful.
        self.assertGreater(outcomes["ok"], 10)
        self.assertGreater(outcomes["error"], 10)


if __name__ == "__main__":
    unittest.main()

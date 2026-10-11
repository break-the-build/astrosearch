"""Signal diagnostics: depth provenance, event-pattern grouping, variability
conflicts. Unit tests use hand-built rows; synthetic end-to-end cases run the
real astrosearch binary with Profile J and five search passes."""
import csv
import math
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import signal_diagnostics as cli  # noqa: E402
from vetting import signal_diagnostics as sd  # noqa: E402

BINARY = Path(os.environ.get("ASTROSEARCH_BIN", ROOT / "target/release/astrosearch"))
csv.field_size_limit(2**31 - 1)


def row(rank, period, epoch, duration, depth, note, events, cls="periodic_dip", score=30.0, odd_even=""):
    evidence = ";".join(f"{e:.8f}:{depth:.4f}:3.0:5" for e in events)
    return {"rank": str(rank), "class": cls, "period_days": str(period), "epoch": str(epoch),
            "duration_days": str(duration), "depth_ppm": str(depth), "amplitude_ppm": str(depth),
            "score": str(score), "odd_even_sigma": odd_even, "event_evidence": evidence, "note": note}


class ProvenanceParsingTests(unittest.TestCase):
    def test_estimators_are_read_from_the_note(self):
        events = [1.0 + 0.5 * k for k in range(10)]
        rows = [
            row(1, 0.5, 1.0, 0.02, 1.9, "x; found in search pass 5 with 55.2% of samples masked by earlier "
                "passes; depth remeasured after transit-masked detrending (box depth 2837.1 ppm)", events),
            row(2, 0.5, 1.001, 0.02, 2694.9, "x; transit-masked local depth -115.7 ppm is not positive; "
                "box depth kept", events),
            row(3, 0.5, 1.0, 0.02, 1500.0, "transit-like", events),
        ]
        a, b, c = sd.from_astrosearch(rows)
        self.assertEqual((a.box_depth_ppm, a.local_depth_ppm, a.search_pass, a.masked_percent),
                         (2837.1, 1.9, 5, 55.2))
        self.assertEqual((b.box_depth_ppm, b.local_depth_ppm), (2694.9, -115.7))
        self.assertTrue(b.depth_estimator.startswith("box"))
        self.assertEqual((c.box_depth_ppm, c.local_depth_ppm, c.depth_estimator), (1500.0, None, "box"))

    def test_wt682_like_rows_form_one_pattern_with_inconsistent_depths(self):
        events = [1.0 + 0.1877 * k for k in range(200)]
        rows = [
            row(21, 0.18772, 1.0, 0.0147, 1.9, "depth remeasured after transit-masked detrending "
                "(box depth 2837.1 ppm)", events),
            row(22, 0.18769, 1.003, 0.0147, 2694.9, "transit-masked local depth -115.7 ppm is not positive; "
                "box depth kept", events[::1][:150]),
            row(24, 0.37544, 1.06, 0.0205, 749.4, "depth remeasured after transit-masked detrending "
                "(box depth 4145.4 ppm)", [1.06 + 0.37544 * k for k in range(100)]),
            row(1, 0.37505, 1.0, 0.0, 0.0, "", [], cls="periodic_variability", score=199.0),
        ]
        report = sd.diagnose(rows)
        groups = report["event_groups"]
        first = next(g for g in groups if any(m["detection"] == "astrosearch:21" for m in g["members"]))
        self.assertEqual({m["detection"] for m in first["members"]}, {"astrosearch:21", "astrosearch:22"})
        self.assertEqual(first["relations"][0]["period_relation"], "same-period")
        check = first["depth_check"]
        self.assertEqual(check["status"], "INCONSISTENT")
        self.assertTrue(any("1.9-2694.9 ppm" in issue for issue in check["issues"]))
        # Every value is still reported.
        self.assertEqual(sorted(m["reported_depth_ppm"] for m in first["members"]), [1.9, 2694.9])
        conflicts = {c["detection"]: c for c in report["variability_conflicts"]}
        self.assertEqual(conflicts["astrosearch:21"]["related_variability"][0]
                         ["transit_period_over_variability_period"], 0.5)
        self.assertEqual(conflicts["astrosearch:24"]["related_variability"][0]
                         ["transit_period_over_variability_period"], 1.0)
        self.assertIn("not tested", conflicts["astrosearch:24"]["status"])

    def test_relations_and_event_matching(self):
        self.assertEqual(sd.period_relation(1.0, 2.0005), "2:1 harmonic")
        self.assertEqual(sd.period_relation(2.0, 3.0), "3:2 harmonic")
        self.assertTrue(sd.period_relation(1.0, 1.37).startswith("unrelated"))
        half = sd.Detection("x", "1", 1.0, 0.0, 0.05, 1.0, "box", 1.0, None, None, None,
                            tuple(float(k) for k in range(20)))
        full = sd.Detection("x", "2", 2.0, 0.0, 0.05, 1.0, "box", 1.0, None, None, None,
                            tuple(float(k) for k in range(0, 20, 2)))
        self.assertEqual(sd.shared_events(half, full), 10)  # every 2P event is a P/2 event


def synthetic(kind: str, seed: int) -> tuple[np.ndarray, np.ndarray]:
    """Two 13-day segments at 10-minute cadence with 400 ppm white noise."""
    rng = np.random.default_rng(seed)
    step = 10 / 1440
    time = np.concatenate([np.arange(1000.0, 1013.0, step), np.arange(1015.0, 1028.0, step)])
    flux = 1.0 + 4e-4 * rng.standard_normal(time.size)

    def box(period, epoch, duration, depth, phase_shift=0.0):
        phase = ((time - epoch - phase_shift * period + period / 2) % period) - period / 2
        return np.where(np.abs(phase) <= duration / 2, depth, 0.0)

    if kind == "transit":
        flux -= box(2.3, 1000.7, 0.09, 2500e-6)
    elif kind == "sinusoid":
        flux += 3000e-6 * np.sin(2 * np.pi * time / 0.6)
    elif kind == "eb":  # unequal primary and secondary eclipses
        flux -= box(1.8, 1000.4, 0.1, 9000e-6) + box(1.8, 1000.4, 0.1, 3000e-6, phase_shift=0.5)
    elif kind == "blend":  # the same binary diluted tenfold by a brighter neighbour
        flux -= 0.1 * (box(1.8, 1000.4, 0.1, 9000e-6) + box(1.8, 1000.4, 0.1, 3000e-6, phase_shift=0.5))
    elif kind != "noise":
        raise ValueError(kind)
    return time, flux


@unittest.skipUnless(BINARY.exists(), f"astrosearch binary not found at {BINARY}")
class SyntheticProfileJTests(unittest.TestCase):
    """End to end: Profile J, five passes, then the diagnostics CLI."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.reports = {}
        for seed, kind in enumerate(("transit", "sinusoid", "eb", "blend", "noise")):
            directory = Path(cls.tmp.name) / kind
            directory.mkdir()
            time, flux = synthetic(kind, 100 + seed)
            curve = directory / "lc.csv"
            np.savetxt(curve, np.column_stack([time, flux, np.full(time.size, 4e-4)]), delimiter=",",
                       header="time,flux,flux_err", comments="", fmt="%.8f")
            candidates = directory / "candidates.csv"
            subprocess.run([str(BINARY), "--input", str(curve), "--output", str(candidates),
                            "--search-profile", "fast", "--transit-search-passes", "5",
                            "--min-period", "0.3", "--max-period", "6", "--threads", "2"],
                           check=True, capture_output=True, text=True)
            output = directory / "diagnostics.json"
            cli.main(["--candidates", str(candidates), "--lightcurve", str(curve), "--output", str(output)])
            import json
            with candidates.open() as stream:
                cls.reports[kind] = (json.loads(output.read_text()), list(csv.DictReader(stream)))

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def best_transit(self, kind, period, tolerance=0.02):
        _, rows = self.reports[kind]
        return next(r for r in rows if r["class"] == "periodic_dip"
                    and abs(float(r["period_days"]) - period) / period < tolerance)

    def member(self, kind, rank):
        report, _ = self.reports[kind]
        for group in report["event_groups"]:
            for m in group["members"]:
                if m["detection"] == f"astrosearch:{rank}":
                    return group, m
        raise AssertionError(rank)

    def test_inputs_are_not_modified(self):
        for kind, (report, rows) in self.reports.items():
            self.assertIn("report-only", report["provenance"], kind)
            with (Path(self.tmp.name) / kind / "candidates.csv").open() as stream:
                self.assertEqual(len(rows), len(list(csv.DictReader(stream))))

    def test_injected_transit_is_consistent_and_not_explained_by_modulation(self):
        transit = self.best_transit("transit", 2.3)
        group, member = self.member("transit", transit["rank"])
        self.assertEqual(member["search_pass"], 1)
        self.assertGreater(member["unmasked_depth_ppm"], 0.6 * float(transit["depth_ppm"]))
        self.assertFalse(any(transit_rank_issue(issue, transit["rank"]) for issue in group["depth_check"]["issues"]),
                         group["depth_check"])
        report, _ = self.reports["transit"]
        conflict = next(c for c in report["variability_conflicts"] if c["detection"] == f"astrosearch:{transit['rank']}")
        verdict = conflict.get("modulation_test", {}).get("verdict", conflict["status"])
        self.assertNotIn("competing explanation", verdict)

    def test_pure_sinusoid_dips_are_explained_by_the_modulation(self):
        report, rows = self.reports["sinusoid"]
        self.assertTrue(any(r["class"] == "periodic_variability" and abs(float(r["period_days"]) - 0.6) < 0.01
                            for r in rows))
        tested = [c for c in report["variability_conflicts"] if c.get("modulation_test", {}).get("status") == "OK"]
        self.assertTrue(tested, report["variability_conflicts"])
        positive = [c for c in tested if (c["modulation_test"]["retained_fraction"] or 0) > 0]
        # Every related dip with a positive unmasked depth loses most of it
        # once the modulation is removed.
        for c in tested:
            retained = c["modulation_test"]["retained_fraction"]
            self.assertTrue(retained is None or retained < sd.EXPLAINED_BELOW, c)
        self.assertGreaterEqual(len(tested), 1)
        self.assertLessEqual(len(positive), len(tested))

    def eclipse_rows(self, kind):
        _, rows = self.reports[kind]
        dips = [r for r in rows if r["class"] == "periodic_dip"]
        primary = next(r for r in dips if abs(float(r["period_days"]) - 1.8) < 0.02)
        secondary = next(r for r in dips if abs(float(r["period_days"]) - 0.9) < 0.01)
        return primary, secondary

    def test_eclipsing_binary_secondary_is_paired_with_the_primary(self):
        # Measured with Profile J, five passes: pass 1 finds the 1.8-day
        # primary; with its eclipses masked, pass 2 reports the secondary as a
        # 0.9-day row whose events all fall on one parity.
        for kind, scale in (("eb", 1.0), ("blend", 0.1)):
            primary, secondary = self.eclipse_rows(kind)
            report, _ = self.reports[kind]
            _, p_member = self.member(kind, primary["rank"])
            _, s_member = self.member(kind, secondary["rank"])
            self.assertEqual(p_member["parity_coverage"], "both parities", kind)
            self.assertTrue(s_member["parity_coverage"].endswith("(effective period 2P)"), s_member)
            self.assertEqual(s_member["search_pass"], 2, s_member)
            pair = next(x for x in report["interleaved_pairs"]
                        if set(x["detections"]) == {f"astrosearch:{primary['rank']}",
                                                    f"astrosearch:{secondary['rank']}"})
            self.assertLess(abs(pair["phase_of_second_events"] - 0.5), sd.SECONDARY_PHASE_WINDOW, pair)
            self.assertIn("eclipsing binary", pair["interpretation"])
            # Injected 9000 and 3000 ppm (x scale); depths stay within 25%.
            depths = sorted(pair["box_depths_ppm"])
            self.assertLess(abs(depths[0] - 3000 * scale) / (3000 * scale), 0.25, pair)
            self.assertLess(abs(depths[1] - 9000 * scale) / (9000 * scale), 0.25, pair)

    def test_eclipsing_binary_dips_are_not_explained_by_modulation(self):
        report, _ = self.reports["eb"]
        tested = [c for c in report["variability_conflicts"] if c.get("modulation_test", {}).get("status") == "OK"]
        for c in tested:  # sharp eclipses are not absorbed by two harmonics
            if c["modulation_test"]["depth_before_ppm"] > 2000:
                self.assertGreater(c["modulation_test"]["retained_fraction"], sd.SURVIVES_ABOVE, c)

    def test_diluted_binary_and_noise_run_cleanly(self):
        for kind in ("blend", "noise"):
            report, rows = self.reports[kind]
            self.assertEqual(report["detections"],
                             sum(r["class"] in ("periodic_dip", "duo_transit") for r in rows), kind)
            for group in report["event_groups"]:
                for m in group["members"]:
                    self.assertIn("unmasked_depth_ppm", m)


def transit_rank_issue(issue: str, rank: str) -> bool:
    return issue.startswith(f"astrosearch:{rank}:")


class ModulationTestUnits(unittest.TestCase):
    def test_box_survives_and_sinusoid_does_not(self):
        rng = np.random.default_rng(3)
        time = np.arange(0, 20, 0.01)
        det = sd.Detection("x", "1", 1.0, 0.25, 0.08, 3000.0, "box", 3000.0, None, None, None, ())
        wave = 1 + 3e-3 * np.cos(2 * np.pi * (time - 0.25)) * -1 + 2e-4 * rng.standard_normal(time.size)
        result = sd.modulation_test(time, wave, det, 1.0)
        self.assertLess(result["retained_fraction"], sd.EXPLAINED_BELOW, result)
        phase = ((time - 0.25 + 0.5) % 1.0) - 0.5
        transit = 1 - 3e-3 * (np.abs(phase) <= 0.04) + 2e-4 * rng.standard_normal(time.size)
        result = sd.modulation_test(time, transit, det, 1.0)
        self.assertGreater(result["retained_fraction"], sd.SURVIVES_ABOVE, result)
        self.assertTrue(math.isclose(sd.unmasked_depth(time, transit, det), 3000, rel_tol=0.1))

    def test_cli_rejects_bad_paths(self):
        with self.assertRaises(SystemExit):
            cli.main(["--candidates", "x.csv", "--output", "x.md"])
        with self.assertRaises(SystemExit):
            cli.main(["--candidates", "x.csv", "--tls", "t.csv", "--output", "r.json"])


if __name__ == "__main__":
    unittest.main()

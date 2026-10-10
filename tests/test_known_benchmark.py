import csv
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from build_known_answer_benchmark import (archive_time_to_btjd, eligible_planets,
                                          select_rows, split_for_host)
from calibrate_detection_threshold import case_values, required_cases
from compare_detection_statistics import paired_summary
from evaluate_benchmark import choose_threshold, close_ephemeris, ephemerides


class KnownAnswerBenchmarkTests(unittest.TestCase):
    def test_archive_planets_are_deduplicated_by_host_and_filtered(self):
        rows = [
            {"pl_name": "Planet b", "hostname": "Star", "tic_id": "123.0",
             "pl_orbper": "2", "pl_orbpererr1": "0.0001",
             "pl_tranmid": "2459000", "pl_tranmiderr1": "0.001", "pl_trandur": "2.4",
             "tran_flag": "1", "pl_controv_flag": "0"},
            {"pl_name": "Planet c", "hostname": "Star", "tic_id": "123",
             "pl_orbper": "4", "pl_orbpererr1": "0.0001",
             "pl_tranmid": "2459001", "pl_tranmiderr1": "0.001", "pl_trandur": "3",
             "tran_flag": "1", "pl_controv_flag": "0"},
            {"pl_name": "Controversial", "tic_id": "456", "pl_orbper": "2",
             "pl_tranmid": "2459000", "pl_trandur": "2", "tran_flag": "1",
             "pl_controv_flag": "1"},
            {"pl_name": "TTV system", "tic_id": "789", "pl_orbper": "2",
             "pl_orbpererr1": "0.0001", "pl_tranmid": "2459000",
             "pl_tranmiderr1": "0.001", "pl_trandur": "2", "tran_flag": "1",
             "pl_controv_flag": "0", "ttv_flag": "1"},
        ]
        selected = eligible_planets(rows)
        self.assertEqual(len(selected), 1)
        self.assertEqual(selected[0]["tic_id"], "123")
        self.assertIn(selected[0]["duration_days"], {0.1, 0.125})

    def test_archive_julian_epoch_converts_to_btjd(self):
        self.assertEqual(archive_time_to_btjd(2457001.25), 1.25)
        self.assertEqual(archive_time_to_btjd(123.25), 123.25)

    def test_split_is_stable_and_sampler_returns_requested_host_count(self):
        self.assertEqual(split_for_host(42, "123"), split_for_host(42, "123"))
        rows = []
        for i in range(100):
            tic = str(10000 + i)
            rows.append({"tic_id": tic, "stratum_period": "short", "pl_name": f"P{i}"})
        selected = select_rows(rows, 20, 91)
        self.assertEqual(len(selected), 40)
        self.assertEqual(sum(row["split"] == "blind" for row in selected), 20)
        self.assertEqual(len({row["tic_id"] for row in selected}), 40)

    def test_threshold_is_selected_only_from_calibration_control_scores(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cases = []
            for i in range(20):
                candidate = root / f"{i}.csv"
                with candidate.open("w", newline="") as stream:
                    writer = csv.writer(stream)
                    writer.writerow(["class", "score"])
                    writer.writerow(["periodic_dip", i + 1])
                cases.append({"case_id": str(i), "truth_class": "inverted_null",
                              "candidates_csv": candidate.name})
            threshold, summary = choose_threshold(cases, root, 0.1, 0)
            self.assertEqual(threshold, 19.0)
            self.assertEqual(summary["triggered_controls"], 2)
            self.assertEqual(summary["control_cases"], 20)

    def test_jackknife_threshold_uses_multi_event_statistic(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cases = []
            rows = [
                (100, "100:1:99:2;101:1:1:2"),
                (40, "100:1:10:2;101:1:10:2"),
            ]
            for i, (score, evidence) in enumerate(rows):
                candidate = root / f"{i}.csv"
                with candidate.open("w", newline="") as stream:
                    writer = csv.writer(stream)
                    writer.writerow(["class", "score", "event_evidence"])
                    writer.writerow(["periodic_dip", score, evidence])
                cases.append({"case_id": str(i), "truth_class": "inverted_null",
                              "candidates_csv": candidate.name})
            threshold, summary = choose_threshold(cases, root, 0.5, 0,
                                                  "jackknife_score")
            self.assertAlmostEqual(threshold, 40 / 2**0.5)
            self.assertEqual(summary["triggered_controls"], 1)

    def test_planet_recovery_uses_best_matching_candidate_statistic(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            candidate = root / "planet.csv"
            with candidate.open("w", newline="") as stream:
                writer = csv.writer(stream)
                writer.writerow(["class", "score", "period_days", "epoch",
                                 "duration_days", "event_evidence"])
                writer.writerow(["periodic_dip", 100, 3, 100, 0.2,
                                 "100:1:99:2;103:1:1:2"])
                writer.writerow(["periodic_dip", 40, 3, 100, 0.2,
                                 "100:1:10:2;103:1:10:2"])
            case = {"truth_class": "confirmed_planet", "baseline_days": "30",
                    "period_days": "3", "epoch_days": "100",
                    "duration_days": "0.2", "candidates_csv": candidate.name}
            _, matched = case_values(case, root, "jackknife_score")
            self.assertAlmostEqual(matched, 40 / 2**0.5)

    def test_required_results_cannot_be_silently_dropped(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name in ("control.csv", "planet.csv", "blind-null.csv"):
                (root / name).write_text("class,score\nperiodic_dip,10\n")
            rows = [
                {"case_id": "cal-pos", "split": "calibration",
                 "truth_class": "confirmed_planet", "candidates_csv": "not-run.csv"},
                {"case_id": "cal-null", "split": "calibration",
                 "truth_class": "inverted_null", "candidates_csv": "control.csv"},
                {"case_id": "blind-pos", "split": "blind",
                 "truth_class": "confirmed_planet", "candidates_csv": "planet.csv"},
                {"case_id": "blind-null", "split": "blind",
                 "truth_class": "inverted_null", "candidates_csv": "blind-null.csv"},
            ]
            required, skipped = required_cases(rows, root)
            self.assertEqual(len(required), 3)
            self.assertEqual(len(skipped), 1)
            (root / "blind-null.csv").unlink()
            with self.assertRaisesRegex(ValueError, "missing results for 1 required"):
                required_cases(rows, root)

    def test_paired_statistics_preserve_case_matching(self):
        result = paired_summary({"a": True, "b": True, "c": False},
                                {"a": False, "b": True, "c": True},
                                bootstrap_replicates=1000, seed=7)
        self.assertEqual(result["first_only"], 1)
        self.assertEqual(result["second_only"], 1)
        self.assertEqual(result["both"], 1)
        self.assertEqual(result["neither"], 0)
        self.assertEqual(result["mcnemar_exact_two_sided_p"], 1.0)

    def test_truth_matching_accepts_preserved_alias_ephemerides(self):
        row = {"period_days": "2.0", "epoch": "100.0", "duration_days": "0.1",
               "alias_ephemerides": "3.000:101.0:0.2:12.0:3:101+104+107"}
        aliases = ephemerides(row)
        self.assertEqual(len(aliases), 2)
        truth = {"period_days": "3.0", "epoch_days": "101", "duration_days": "0.2"}
        self.assertTrue(close_ephemeris(aliases[1], truth, 30, 1.0))


if __name__ == "__main__":
    unittest.main()

import csv
import json
import sqlite3
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
import sys
sys.path.insert(0, str(ROOT / "scripts"))

from candidate_registry import ingest, transition


class CandidateRegistryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.database = self.root / "candidates.sqlite"
        self.lightcurve = self.root / "curve.csv"
        self.lightcurve.write_text("time,flux\n1,1\n2,1\n")
        self.candidates = self.root / "candidates.csv"
        with self.candidates.open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["rank", "class", "score", "period_days", "epoch",
                             "duration_days", "depth_ppm", "event_evidence"])
            writer.writerow([1, "periodic_dip", 35.0, 3.2, 100.0, 0.1, 500,
                             "100:500:12:3;103.2:500:11:3"])

    def tearDown(self):
        self.temp.cleanup()

    def test_ingest_is_idempotent_and_starts_at_detected_signal(self):
        first = ingest(self.database, "TIC123", self.lightcurve,
                       self.candidates, "abc123")
        repeat = ingest(self.database, "TIC123", self.lightcurve,
                        self.candidates, "abc123")
        self.assertTrue(first["inserted"])
        self.assertFalse(repeat["inserted"])
        with sqlite3.connect(self.database) as connection:
            self.assertEqual(connection.execute("SELECT COUNT(*) FROM runs").fetchone()[0], 1)
            detection = connection.execute("SELECT detection_id,status FROM detections").fetchone()
            self.assertEqual(detection[1], "DETECTED_SIGNAL")
        self.detection_id = detection[0]

    def test_human_evidence_and_sequential_promotion_are_required(self):
        result = ingest(self.database, "TIC123", self.lightcurve,
                        self.candidates, "abc123")
        with self.assertRaisesRegex(ValueError, "actor and evidence"):
            transition(self.database, result["run_id"], "VETTED_CANDIDATE", "", "")
        with sqlite3.connect(self.database) as connection:
            detection_id = connection.execute("SELECT detection_id FROM detections").fetchone()[0]
        with self.assertRaisesRegex(ValueError, "advance one stage"):
            transition(self.database, detection_id, "CONFIRMED_DISCOVERY",
                       "reviewer", "paper and follow-up")
        promoted = transition(self.database, detection_id, "VETTED_CANDIDATE",
                              "reviewer", "vetting/report-001")
        self.assertEqual(promoted["previous_status"], "DETECTED_SIGNAL")
        with sqlite3.connect(self.database) as connection:
            status = connection.execute("SELECT status FROM detections WHERE detection_id=?",
                                        (detection_id,)).fetchone()[0]
            history = connection.execute("SELECT COUNT(*) FROM status_history WHERE detection_id=?",
                                         (detection_id,)).fetchone()[0]
        self.assertEqual(status, "VETTED_CANDIDATE")
        self.assertEqual(history, 2)

    def test_ingest_preserves_manifest_with_sector_product_metadata(self):
        manifest = self.root / "manifest.json"
        manifest.write_text(json.dumps({"tic_id": 123, "products": [
            {"sector": 94, "product_uri": "mast:hlsp/qlp/s0094"}]}))
        result = ingest(self.database, "TIC123", self.lightcurve,
                        self.candidates, "abc123", manifest)
        with sqlite3.connect(self.database) as connection:
            stored = connection.execute(
                "SELECT provenance_json FROM runs WHERE run_id=?", (result["run_id"],)
            ).fetchone()[0]
        provenance = json.loads(stored)
        self.assertEqual(provenance["manifest"]["products"][0]["sector"], 94)
        self.assertTrue(provenance["manifest_sha256"])


if __name__ == "__main__":
    unittest.main()

import sys
import tempfile
import unittest
from types import SimpleNamespace
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from candidate_preflight import select_tce_files
from survey_ledger import (claim_target, connect, begin_run,
                           config_digest, finish_run, record_target,
                           should_defer)
from tess_survey import RequestTimeoutAdapter, atomic_write_text, retry_network
import requests


class SreOperationalTests(unittest.TestCase):
    def test_config_fingerprint_changes_with_candidate_grouping(self):
        config = SimpleNamespace(
            authors=["QLP", "SPOC"], min_sector=90, min_exptime=200,
            min_baseline=1095, min_teff=2700, max_teff=5200,
            max_distance_pc=100, max_radius=.9, max_tmag=14,
            max_sector_noise_ratio=5.0, min_period=.15, max_period=30,
            periods=None, threads=0, trend_method="median",
            candidate_grouping="ephemeris")
        original = config_digest(config)
        config.candidate_grouping = "events"
        self.assertNotEqual(original, config_digest(config))

    def test_tce_range_is_selected_when_requested_sector_is_inside_range(self):
        files = select_tce_files(
            ["/tess/tess-s0001-s0055_dvr-tcestats.csv",
             "/tess/tess-s0056-s0056_dvr-tcestats.csv"],
            [32], "https://archive.stsci.edu/tess/bulk_downloads/index.html")
        self.assertEqual(len(files), 1)
        self.assertEqual(next(iter(files.values())), (1, 55))

    def test_atomic_claim_prevents_two_runs_from_owning_one_target(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "ledger.sqlite"
            first, second = connect(path), connect(path)
            begin_run(first, "run-1", path, "local", "config")
            begin_run(second, "run-2", path, "local", "config")
            row = {"tic_id": "123", "teff_k": "3500"}
            self.assertTrue(claim_target(first, "run-1", "123", row, "config"))
            self.assertFalse(claim_target(second, "run-2", "123", row, "config"))
            first.close(); second.close()

    def test_expired_target_claim_can_be_recovered(self):
        with tempfile.TemporaryDirectory() as directory:
            db = connect(Path(directory) / "ledger.sqlite")
            row = {"tic_id": "123", "teff_k": "3500"}
            begin_run(db, "run-1", Path(directory), "local", "config")
            begin_run(db, "run-2", Path(directory), "local", "config")
            self.assertTrue(claim_target(db, "run-1", "123", row, "config"))
            db.execute("UPDATE target_state SET updated_utc='2000-01-01T00:00:00+00:00' WHERE tic_id='123'")
            db.commit()
            self.assertTrue(claim_target(db, "run-2", "123", row, "config"))
            db.close()

    def test_interrupted_run_releases_target_for_immediate_retry(self):
        with tempfile.TemporaryDirectory() as directory:
            db = connect(Path(directory) / "ledger.sqlite")
            row = {"tic_id": "123", "teff_k": "3500"}
            begin_run(db, "run-1", Path(directory), "local", "config")
            self.assertTrue(claim_target(db, "run-1", "123", row, "config"))
            finish_run(db, "run-1", "interrupted")
            saved = db.execute("SELECT status FROM target_state WHERE tic_id='123'").fetchone()
            self.assertEqual(saved["status"], "error")
            self.assertFalse(should_defer(db, "123", row, "config", 30, False))
            db.close()

    def test_completed_no_data_is_deferred_but_errors_are_retryable(self):
        with tempfile.TemporaryDirectory() as directory:
            db = connect(Path(directory) / "ledger.sqlite")
            row = {"tic_id": "123", "teff_k": "3500"}
            begin_run(db, "run-1", Path(directory), "local", "config")
            record_target(db, "run-1", "123", row, "config", {"status": "completed_no_data"})
            self.assertTrue(should_defer(db, "123", row, "config", 30, False))
            record_target(db, "run-1", "123", row, "config", {"status": "error"})
            self.assertFalse(should_defer(db, "123", row, "config", 30, False))
            db.close()

    def test_retry_network_retries_connection_failure_only_with_bounded_backoff(self):
        attempts = 0
        def operation():
            nonlocal attempts
            attempts += 1
            if attempts == 1:
                raise requests.exceptions.ConnectionError("temporary")
            return "ok"
        with patch("tess_survey.time.sleep"):
            self.assertEqual(retry_network(operation, "test"), "ok")
        self.assertEqual(attempts, 2)

    def test_mast_adapter_injects_timeout_but_preserves_explicit_override(self):
        adapter = RequestTimeoutAdapter(12)
        request = requests.Request("GET", "https://mast.stsci.edu").prepare()
        with patch("requests.adapters.HTTPAdapter.send", return_value=None) as send:
            adapter.send(request)
            self.assertEqual(send.call_args.kwargs["timeout"], 12)
            adapter.send(request, timeout=3)
            self.assertEqual(send.call_args.kwargs["timeout"], 3)

    def test_atomic_text_replaces_complete_file_and_leaves_no_temp(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "manifest.json"
            output.write_text("old")
            atomic_write_text(output, "new")
            self.assertEqual(output.read_text(), "new")
            self.assertEqual([p.name for p in Path(directory).iterdir()], [output.name])


if __name__ == "__main__":
    unittest.main()

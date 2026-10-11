"""TLS iteration bookkeeping: observed transit times, event-based relations
between iterations, and alias labels. Uses TLS-shaped result objects, so it
runs without the optional transitleastsquares dependency."""
import sys
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import tls_search  # noqa: E402


def row(iteration, period, epoch, duration, times):
    return {"iteration": iteration, "period_days": period, "epoch_days": epoch, "duration_days": duration,
            "depth_ppm": 3000.0, "observed_transit_times": ";".join(f"{t:.8f}" for t in times)}


class TlsBookkeepingTests(unittest.TestCase):
    def test_observed_times_keep_transits_with_data(self):
        result = SimpleNamespace(transit_times=np.array([1.0, 2.0, 3.0, 4.0]),
                                 per_transit_count=np.array([5, 0, 3, 0]))
        self.assertEqual(tls_search.observed_times(result), [1.0, 3.0])
        self.assertEqual(tls_search.observed_times(SimpleNamespace()), [])

    def test_alias_labels(self):
        prior = [{"period_days": 0.37544}]
        self.assertEqual(tls_search.alias_label(0.375441, prior), "same-period")
        self.assertEqual(tls_search.alias_label(0.75088, prior), "integer")
        self.assertEqual(tls_search.alias_label(0.18772, prior), "integer")
        self.assertEqual(tls_search.alias_label(0.56316, prior), "half-integer")
        self.assertEqual(tls_search.alias_label(0.4, prior), "none")

    def test_repeated_feature_is_reported_by_shared_events(self):
        # WT 682 pattern: the same 0.37544-day feature at a slightly shifted epoch.
        first = row(1, 0.3754436, 1654.150625, 0.0019, [1654.150625 + 0.3754436 * k for k in range(300)])
        again = row(2, 0.3754407, 1654.157392, 0.0019, [1654.157392 + 0.3754407 * k for k in range(290)])
        relation = tls_search.event_relations(again, [first])
        self.assertTrue(relation.startswith("1:"), relation)
        shared, total = relation.split()[0].split(":")[1].split("/")
        self.assertGreater(int(shared) / int(total), 0.5, relation)
        self.assertIn("same-period", relation)
        unrelated = row(3, 2.17, 1650.3, 0.08, [1650.3 + 2.17 * k for k in range(40)])
        self.assertEqual(tls_search.event_relations(unrelated, [first]), "")

    def test_depth_is_converted_from_bottom_flux(self):
        # WT 682 iteration 1: TLS depth 0.99613537 is a 3,864.6 ppm transit.
        self.assertAlmostEqual(tls_search.depth_ppm_from_tls(0.99613537), 3864.63, places=2)

    def test_mask_width_is_unchanged(self):
        # Search results must match the previous release; see MASK_DURATIONS.
        self.assertEqual(tls_search.MASK_DURATIONS, 1.0)


if __name__ == "__main__":
    unittest.main()

"""Offline tests of the TESS export sidecars the gate relies on: quality-bit
decoding (momentum dumps, scattered light) and the crowding-correction flag
that decides whether vetting applies the TIC dilution factor."""
import sys
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

try:
    import export_tess  # needs lightkurve (locked environment)
except ImportError:  # pragma: no cover
    export_tess = None


def fake_curve(time, quality, meta=None):
    return SimpleNamespace(time=SimpleNamespace(value=np.asarray(time, float)),
                           quality=np.asarray(quality, int), meta=meta or {})


@unittest.skipIf(export_tess is None, "lightkurve not installed")
class ExportMetadataTests(unittest.TestCase):
    def test_spoc_quality_bits_decode_to_dumps_and_scattered_light_spans(self):
        time = 1400 + np.arange(20) * 0.01
        quality = np.zeros(20, int)
        quality[3] = 32                    # momentum dump
        quality[10:13] = 4096              # scattered light, contiguous
        quality[15] = 2048 | 32            # both
        quality[18] = 1024                 # unrelated bit
        events = export_tess.quality_events(fake_curve(time, quality), "SPOC")
        np.testing.assert_allclose(events["momentum_dump_times"], [time[3], time[15]])
        spans = events["scattered_light_intervals"]
        self.assertEqual(len(spans), 2)
        np.testing.assert_allclose(spans[0], [time[10], time[12]])
        np.testing.assert_allclose(spans[1], [time[15], time[15]])

    def test_non_spoc_quality_is_not_interpreted(self):
        events = export_tess.quality_events(fake_curve([1, 2], [32, 4096]), "QLP")
        self.assertNotIn("momentum_dump_times", events)

    def test_crowding_flag_only_for_spoc_pdcsap_with_crowdsap(self):
        rows = np.array([[1400.0, 1.0, 1e-3], [1401.0, 1.0, 1e-3]])
        raw = fake_curve([1400.0, 1401.0], [0, 0])

        def flag(author, meta):
            lc = fake_curve([1400.0, 1401.0], [0, 0], meta)
            return export_tess.metadata(lc, raw, "TIC 1", 1, author, 120, rows)["crowding_corrected"]

        spoc = {"FLUX_ORIGIN": "pdcsap_flux", "CROWDSAP": 0.98}
        self.assertTrue(flag("SPOC", spoc))
        self.assertTrue(flag("TESS-SPOC", spoc))
        self.assertIsNone(flag("SPOC", {"FLUX_ORIGIN": "sap_flux", "CROWDSAP": 0.98}))
        self.assertIsNone(flag("SPOC", {"FLUX_ORIGIN": "pdcsap_flux"}))
        self.assertIsNone(flag("QLP", {"FLUX_ORIGIN": "kspsap_flux", "CROWDSAP": 0.98}))


if __name__ == "__main__":
    unittest.main()

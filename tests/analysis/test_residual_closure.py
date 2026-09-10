"""Known monotone and periodic signals distinguish drift from signed cycling."""
import sys
import unittest
from pathlib import Path
try:
    import numpy as np
except ModuleNotFoundError:
    np = None
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
if np is not None:
    from analyze_residual_closure import temporal, compact_balance


@unittest.skipIf(np is None, "C2b postprocessing tests require the project analysis .venv")
class ResidualEvidenceTests(unittest.TestCase):
    def test_monotone_decay_has_no_sign_cycle(self):
        t = np.linspace(0, 10, 1001)
        s = temporal(t, 2 - .1*t)
        self.assertEqual(s['sign_changes'], 0)
        self.assertAlmostEqual(s['linear_slope_per_normalized_time'], -.1)
        self.assertLess(s['detrended_stddev'], 1e-14)
        self.assertIsNone(s['sampled_peak_period'])
        self.assertLess(s['last_quarter_rms'], s['first_quarter_rms'])

    def test_signed_periodic_signal_is_resolved(self):
        t = np.arange(1000)*.01
        s = temporal(t, np.cos(2*np.pi*t))
        self.assertEqual(s['sign_changes'], 20)
        self.assertAlmostEqual(s['sampled_peak_period'], 1)
        self.assertGreater(s['sampled_peak_power_fraction'], .99)

    def test_numerical_balance_preserves_distinct_cell_estimate(self):
        source = {'interior_cell_estimate':{'mass_conservation_error': .01},
                  'boundary_fluxes': {'inlet':{'outward_total':[-2,0,0,0]},
                                      'outlet':{'outward_total':[2.000002,0,0,0]}}}
        result = compact_balance(source)
        self.assertEqual(result['interior_cell_estimate']['mass_conservation_error'], .01)
        self.assertAlmostEqual(result['numerical_inlet_exit_mass_mismatch'], 1e-6, places=11)
        self.assertNotIn('numerical_inlet_exit_mass_mismatch', source)


if __name__ == '__main__':
    unittest.main()

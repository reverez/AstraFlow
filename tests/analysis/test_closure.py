import json
import math
from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from analyze_scientific_closure import richardson, stability, norms, read_csv, read_json, interpolate, validate_run
class AnalysisTests(unittest.TestCase):
    def test_known_orders_and_signs(self):
        for p in (1,2,3):
            for sign in (-1,1):
                exact=7
                values=[exact+sign*h**p for h in (0.5,0.25,0.125)]
                r=richardson(*values)
                self.assertTrue(r['valid'])
                self.assertAlmostEqual(r['order'],p)
                self.assertAlmostEqual(r['extrapolated'],exact)
                self.assertAlmostEqual(r['gci_percent'],100*1.25*abs((values[2]-exact)/values[2]))
    def test_reject_invalid_sequences(self):
        for seq in ((1,2,1.5),(1,1,1),(1,2,4),(1,0.5,0)):
            self.assertFalse(richardson(*seq)['valid'])
        with self.assertRaises(ValueError): richardson(1,2,3,r=1)
        with self.assertRaises(ValueError): richardson(1,2,float('nan'))
    def test_window_statistics(self):
        r=stability([0,10,20],[1,2,3])
        self.assertAlmostEqual(r['relative_range'],1)
        self.assertAlmostEqual(r['relative_stddev'],math.sqrt(2/3)/2)
        self.assertAlmostEqual(r['relative_drift'],1)
        self.assertEqual(stability([0,1],[0,0])['relative_range'],0)
        with self.assertRaises(ValueError): stability([1,1],[2,2])
    def test_norms_and_interpolation(self):
        self.assertEqual(norms([1,2,3],[1,2,3]),dict(L1=0,L2=0,Linf=0))
        self.assertEqual(interpolate([0,1],[1,3],[0,0.5,1]),[1,2,3])
        with self.assertRaises(ValueError): interpolate([0,1],[1,3],[-0.1])
    def test_parsing_rejects_incomplete_and_nonfinite(self):
        with tempfile.TemporaryDirectory() as name:
            p=Path(name)
            (p/'history.csv').write_text('iteration,value\n1,2\n2,3\n')
            self.assertEqual(len(read_csv(p/'history.csv')),2)
            (p/'history.csv').write_text('iteration,value\n2,2\n1,3\n')
            with self.assertRaises(ValueError): read_csv(p/'history.csv')
            (p/'x.json').write_text('{"x": NaN}')
            with self.assertRaises(ValueError): read_json(p/'x.json')
            (p/'summary.json').write_text(json.dumps({'termination_reason':'iteration_limit','converged':False}))
            (p/'config.json').write_text('{}')
            for f in ('engineering.csv','convergence.csv'):
                (p/f).write_text('iteration,value\n1,2\n')
            with self.assertRaisesRegex(ValueError,'not steady'): validate_run(p)
if __name__=='__main__': unittest.main()

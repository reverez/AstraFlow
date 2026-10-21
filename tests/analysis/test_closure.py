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
    def test_independent_closure_acceptance_and_tampering(self):
        from analyze_scientific_closure import OBSERVABLES, RESIDUALS
        import csv
        with tempfile.TemporaryDirectory() as name:
            root=Path(name)
            config={'mesh':{'nx':8,'nr':4},'runtime':{'backend':'cpu','precision':'double','residual_tolerance':1e-6},
                    'convergence':{'minimum_iterations':2000,'window_iterations':2000,'sampling_interval':20,
                                   'relative_residual_target':0,'observable_tolerance':5e-4,
                                   'mass_mismatch_tolerance':1e-3,'station_spread_tolerance':2e-3}}
            summary={'grid':[8,4],'backend':'cpu','precision':'double','termination_reason':'steady_converged',
                     'converged':True,'iterations':4000,'final_residuals':[1e-7]*4,'initial_residuals':[.01]*4,
                     **dict.fromkeys(OBSERVABLES,1.0),'exit_axial_velocity':1.,'exit_temperature':1.,
                     'mass_conservation_error':.0001,'mass_flow_spread':.0002}
            (root/'config.json').write_text(json.dumps(config))
            (root/'summary.json').write_text(json.dumps(summary))
            (root/'final_state.vts').write_text('field checked separately by fields()')
            with (root/'engineering.csv').open('w') as f:
                writer=csv.writer(f);writer.writerow(('iteration',*OBSERVABLES,'mass_conservation_error','mass_flow_spread'))
                writer.writerows((i,*([1.0]*6),.0001,.0002) for i in range(2000,4001,20))
            with (root/'convergence.csv').open('w') as f:
                writer=csv.writer(f);writer.writerow(('iteration',*RESIDUALS))
                writer.writerows((i,*([1e-7]*4)) for i in range(2000,4001,20))
            result,st=validate_run(root)
            self.assertTrue(result['converged'])
            self.assertEqual(st['mass_flow']['relative_range'],0)
            summary['mass_flow']=2
            (root/'summary.json').write_text(json.dumps(summary))
            with self.assertRaisesRegex(ValueError,'engineering mismatch'):validate_run(root)
            summary['mass_flow']=1
            (root/'summary.json').write_text(json.dumps(summary))
            lines=(root/'engineering.csv').read_text().splitlines()
            (root/'engineering.csv').write_text('\n'.join(lines[:2]+lines[3:])+'\n')
            with self.assertRaisesRegex(ValueError,'Missing engineering'):validate_run(root)
if __name__=='__main__': unittest.main()

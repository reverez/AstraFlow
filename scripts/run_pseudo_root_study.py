#!/usr/bin/env python3
"""V1.2 secondary study: do physical-time and local pseudo-time converge to the same discrete root?

For each pre-registered residual tolerance the same 128x32 case is integrated twice, once with
SSP-RK2 physical time and once with local pseudo-time RK2, each stopping at the unchanged
convergence monitor. The analysis reports the nondimensional conservative-field Linf difference
between the two endpoints per tolerance; H2 is read from how that difference scales with the
tolerance, never from a single fixed field bound. Completed runs are never overwritten.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from analyze_scientific_closure import read_json, fields, norms, validate_run

ROOT = Path(__file__).resolve().parents[1]
TOLERANCES = (1e-6, 3e-7, 1e-7)
CAP = 400000


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run_case(cli, base, tolerance, pseudo, out):
    config = json.loads(json.dumps(base))
    config['runtime'].update(residual_tolerance=tolerance, max_iterations=CAP)
    config['steady_solver'] = {'enabled': pseudo, 'local_time_step': True, 'cfl': 0.4}
    config['output'].pop('directory', None)
    request = out.with_suffix('.request.json')
    if out.exists() or request.exists():
        print(f'Preserving existing {out}', flush=True)
        return
    request.write_text(json.dumps(config, indent=2) + '\n')
    print(f'Running tolerance={tolerance} {"pseudo" if pseudo else "physical"}: {out}', flush=True)
    with out.with_suffix('.log').open('x') as log:
        code = subprocess.run([str(cli), '--config', str(request), '--output', str(out)],
                              stdout=log, stderr=subprocess.STDOUT).returncode
    meta = {'exit_code': code, 'cli_sha256': digest(cli), 'tolerance': tolerance, 'pseudo': pseudo}
    (out / 'provenance.json').write_text(json.dumps(meta, indent=2) + '\n')


def analyze(root):
    report = {'tolerances': [], 'hypothesis': 'H2: same discrete root',
              'reading': 'pre-registered in docs/verification/v12_closure_protocol.md'}
    linf = {}
    for tolerance in TOLERANCES:
        pair = {}
        for mode in ('physical', 'pseudo'):
            d = root / f'tolerance-{tolerance:g}' / mode
            if not (d / 'summary.json').exists():
                pair[mode] = {'status': 'missing'}
                continue
            s = read_json(d / 'summary.json')
            entry = {'termination': s['termination_reason'], 'iterations': s['iterations'],
                     'wall_ms': read_json(d / 'performance.json')['wall_ms'],
                     'final_residuals': s['final_residuals'],
                     'numerical_boundary_mass_flux_mismatch': s['numerical_boundary_mass_flux_mismatch']}
            try:
                validate_run(d)
                entry['independently_validated'] = True
            except ValueError as error:
                entry['independently_validated'] = False
                entry['validation_error'] = str(error)
            pair[mode] = entry
        both = all(pair[m].get('independently_validated') for m in ('physical', 'pseudo'))
        if both:
            ref = read_json(root / f'tolerance-{tolerance:g}' / 'physical' / 'performance.json')['reference_scales']
            factors = [ref['density'], ref['density'] * ref['velocity'], ref['density'] * ref['velocity'], ref['pressure']]
            _, a = fields(root / f'tolerance-{tolerance:g}' / 'physical')
            _, b = fields(root / f'tolerance-{tolerance:g}' / 'pseudo')
            component = {name: norms([x / f for x in aa], [x / f for x in bb])
                         for name, aa, bb, f in zip(('rho', 'rho_u', 'rho_v', 'rho_E'), a, b, factors)}
            pair['field_difference'] = component
            pair['Linf_all'] = max(v['Linf'] for v in component.values())
            linf[tolerance] = pair['Linf_all']
            sp = read_json(root / f'tolerance-{tolerance:g}' / 'physical' / 'summary.json')
            sq = read_json(root / f'tolerance-{tolerance:g}' / 'pseudo' / 'summary.json')
            pair['engineering_relative_differences'] = {
                k: abs((sp[k] - sq[k]) / sq[k]) for k in
                ('mass_flow', 'estimated_thrust', 'specific_impulse', 'throat_mach', 'exit_mach', 'exit_pressure')}
            pair['wall_time_ratio_physical_over_pseudo'] = pair['physical']['wall_ms'] / pair['pseudo']['wall_ms']
        pair['tolerance'] = tolerance
        report['tolerances'].append(pair)
    ordered = [linf[t] for t in TOLERANCES if t in linf]
    report['complete'] = len(ordered) == len(TOLERANCES)
    if report['complete']:
        monotone = all(x > y for x, y in zip(ordered, ordered[1:]))
        ratio = ordered[0] / ordered[-1] if ordered[-1] > 0 else float('inf')
        report['linf_monotone_decrease'] = monotone
        report['linf_ratio_loosest_over_tightest'] = ratio
        report['H2_supported'] = monotone and ratio >= 3
        report['H2_refuted'] = not monotone
    else:
        report['H2_supported'] = None
        report['inconclusive_reason'] = 'at least one tolerance lacks two independently validated endpoints'
    (root / 'analysis.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=('run', 'analyze'))
    p.add_argument('--root', type=Path, default=ROOT / 'runs/v12-pseudo-root')
    p.add_argument('--config', type=Path, default=ROOT / 'examples/rocket_nozzle/v12_closure.json')
    p.add_argument('--cli', type=Path, default=ROOT / 'build/astraflow_cli')
    args = p.parse_args()
    root = args.root.resolve()
    if args.action == 'run':
        root.mkdir(parents=True, exist_ok=True)
        base = read_json(args.config)
        for tolerance in TOLERANCES:
            (root / f'tolerance-{tolerance:g}').mkdir(exist_ok=True)
            for pseudo in (False, True):
                run_case(args.cli.resolve(), base, tolerance, pseudo,
                         root / f'tolerance-{tolerance:g}' / ('pseudo' if pseudo else 'physical'))
    print(json.dumps(analyze(root), indent=2))
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (ValueError, RuntimeError, OSError) as error:
        print(f'Pseudo root study stopped: {error}', file=sys.stderr)
        sys.exit(2)

#!/usr/bin/env python3
"""Nozzle expansion-ratio sweep: steady 2D axisymmetric thrust against quasi-1D isentropic theory.

Each case changes only the exit radius of the inviscid isentropic-nozzle example and is integrated
to the steady convergence monitor with local pseudo-time. The quasi-1D reference uses the choked
isentropic mass flow and the supersonic area-Mach relation at the exit, with the same ambient
pressure in the thrust equation. Agreement is a consistency check between the 2D solver and the
one-dimensional idealisation; it is not experimental validation. Completed cases are reused.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import subprocess
import sys
from analyze_scientific_closure import read_json

ROOT = Path(__file__).resolve().parents[1]


def area_mach(ratio, gamma, supersonic):
    lo, hi = (1.0, 20.0) if supersonic else (1e-9, 1.0)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        a = (2 / (gamma + 1) * (1 + (gamma - 1) / 2 * mid * mid)) ** ((gamma + 1) / (2 * (gamma - 1))) / mid
        if (a < ratio) == supersonic:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def quasi_1d(config, exit_radius):
    gas, bc, geo = config['gas'], config['boundary_conditions'], config['geometry']
    g, R, p0, t0, pa = gas['gamma'], gas['gas_constant'], bc['P0'], bc['T0'], bc['back_pressure']
    rt, re = geo['throat_radius'], exit_radius
    a_star, a_e = math.pi * rt * rt, math.pi * re * re
    mdot = a_star * p0 / math.sqrt(t0) * math.sqrt(g / R) * ((g + 1) / 2) ** (-(g + 1) / (2 * (g - 1)))
    me = area_mach(a_e / a_star, g, True)
    factor = 1 + (g - 1) / 2 * me * me
    te, pe = t0 / factor, p0 / factor ** (g / (g - 1))
    ue = me * math.sqrt(g * R * te)
    thrust = mdot * ue + (pe - pa) * a_e
    return {'mass_flow': mdot, 'exit_mach': me, 'exit_pressure': pe, 'exit_velocity': ue,
            'thrust': thrust, 'specific_impulse': thrust / (mdot * 9.80665), 'area_ratio': a_e / a_star}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, default=ROOT / 'runs/nozzle-sweep')
    p.add_argument('--config', type=Path, default=ROOT / 'examples/isentropic_nozzle/config.json')
    p.add_argument('--cli', type=Path, default=ROOT / 'build/astraflow_cli')
    p.add_argument('--exit-radii', type=float, nargs='+',
                   default=[0.0055, 0.006, 0.00675, 0.0075, 0.0085, 0.0095])
    p.add_argument('--max-iterations', type=int, default=100000)
    args = p.parse_args()
    root = args.root.resolve()
    root.mkdir(parents=True, exist_ok=True)
    base = read_json(args.config)
    rows = []
    for re in args.exit_radii:
        name = f'exit-{re:.5f}'
        out = root / name
        config = json.loads(json.dumps(base))
        config['geometry']['exit_radius'] = re
        config['numerics']['limiter'] = 'van_albada'
        config['steady_solver'] = {'enabled': True, 'local_time_step': True, 'cfl': 0.4}
        config['runtime'].update(max_iterations=args.max_iterations, residual_tolerance=1e-6)
        config['convergence'] = {'minimum_iterations': 5000, 'sampling_interval': 20, 'window_iterations': 2000,
                                 'relative_residual_target': 0, 'observable_tolerance': 5e-4,
                                 'mass_gate': 'numerical_boundary', 'mass_mismatch_tolerance': 1e-6,
                                 'station_spread_tolerance': 2e-3, 'balance_interval': 0}
        config['output'].pop('directory', None)
        request = root / f'{name}.json'
        if not (out / 'summary.json').exists():
            if out.exists():
                raise RuntimeError(f'Incomplete case directory exists, remove or rename: {out}')
            request.write_text(json.dumps(config, indent=2) + '\n')
            print(f'Running {name}', flush=True)
            with (root / f'{name}.log').open('w') as log:
                subprocess.run([str(args.cli), '--config', str(request), '--output', str(out)],
                               stdout=log, stderr=subprocess.STDOUT, check=False)
        if not (out / 'summary.json').exists():
            rows.append({'exit_radius': re, 'status': 'failed'})
            continue
        s = read_json(out / 'summary.json')
        ref = quasi_1d(config, re)
        row = {'exit_radius': re, 'area_ratio': ref['area_ratio'], 'termination': s['termination_reason'],
               'iterations': s['iterations'], 'converged': s['converged'],
               'thrust_2d': s['estimated_thrust'], 'thrust_1d': ref['thrust'],
               'mass_flow_2d': s['mass_flow'], 'mass_flow_1d': ref['mass_flow'],
               'exit_mach_2d': s['exit_mach'], 'exit_mach_1d': ref['exit_mach'],
               'exit_pressure_2d': s['exit_pressure'], 'exit_pressure_1d': ref['exit_pressure'],
               'isp_2d': s['specific_impulse'], 'isp_1d': ref['specific_impulse'],
               'numerical_boundary_mass_flux_mismatch': s['numerical_boundary_mass_flux_mismatch']}
        for k in ('thrust', 'mass_flow', 'exit_mach', 'exit_pressure', 'isp'):
            row[f'{k}_relative_difference'] = (row[f'{k}_2d'] - row[f'{k}_1d']) / row[f'{k}_1d']
        rows.append(row)
    with (root / 'sweep.csv').open('w') as f:
        keys = sorted({k for r in rows for k in r})
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    converged = [r for r in rows if r.get('converged')]
    report = {'cases': rows, 'converged_cases': len(converged), 'total_cases': len(rows),
              'reference': 'quasi-1D isentropic choked flow with supersonic exit branch; same ambient pressure',
              'max_abs_thrust_relative_difference_converged': max((abs(r['thrust_relative_difference']) for r in converged), default=None),
              'max_abs_mass_flow_relative_difference_converged': max((abs(r['mass_flow_relative_difference']) for r in converged), default=None)}
    (root / 'sweep.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    print(json.dumps(report, indent=2))
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (ValueError, RuntimeError, OSError) as error:
        print(f'Sweep stopped: {error}', file=sys.stderr)
        sys.exit(2)

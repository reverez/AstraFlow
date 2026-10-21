#!/usr/bin/env python3
"""V1.3: continue the V1.2 medium trajectory and run the fine grid with pre-registered budgets.

See docs/verification/v13_continuation_protocol.md. The V1.2 directories are read only; every
V1.3 output has its own directory with provenance. The manifest it writes is consumed by
analyze_scientific_closure.py like a normal closure study.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from analyze_scientific_closure import read_json, fields, validate_run

ROOT = Path(__file__).resolve().parents[1]
MEDIUM_CONTINUATION = 700000
FINE_CAP = 2500000


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def state_file(run, out):
    """Schema-1 physical conservative state rebuilt from a preserved final_state.vts."""
    config = read_json(run / 'config.json')
    summary = read_json(run / 'summary.json')
    _, conservative = fields(run)
    rows = list(zip(*conservative))
    payload = {'schema': 1, 'source_run': str(run), 'source_sha256': {k: digest(run / k) for k in ('config.json', 'summary.json', 'final_state.vts')},
               'iterations': summary['iterations'], 'physical_simulated_time': summary['simulated_time'],
               'configuration': {k: config[k] for k in ('problem', 'gas', 'geometry', 'mesh', 'boundary_conditions', 'numerics')},
               'physical_conservative': [list(r) for r in rows]}
    out.write_text(json.dumps(payload) + '\n')
    return payload


def run_cli(cli, config, out, log, initial=None):
    args = [str(cli), '--config', str(config), '--output', str(out)]
    if initial:
        args += ['--initial-state', str(initial)]
    with log.open('x') as stream:
        return subprocess.run(args, stdout=stream, stderr=subprocess.STDOUT).returncode


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--v12', type=Path, default=ROOT / 'runs/v12-closure')
    p.add_argument('--root', type=Path, default=ROOT / 'runs/v13-closure')
    p.add_argument('--cli', type=Path, default=ROOT / 'build/astraflow_cli')
    p.add_argument('--skip-fine', action='store_true')
    args = p.parse_args()
    v12, root, cli = args.v12.resolve(), args.root.resolve(), args.cli.resolve()
    root.mkdir(parents=True, exist_ok=True)
    manifest = {'schema': 1, 'backend': 'cuda', 'precision': 'double', 'study': 'v13', 'runs': []}
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()

    coarse = v12 / '128x32' / 'attempt-001'
    validate_run(coarse)
    manifest['runs'].append({'grid': [128, 32], 'directory': str(coarse), 'converged': True, 'termination_reason': 'steady_converged', 'source': 'v12'})

    medium_src = v12 / '256x64' / 'attempt-001'
    src_summary = read_json(medium_src / 'summary.json')
    medium_dir = root / '256x64'
    medium_dir.mkdir(exist_ok=True)
    if src_summary['converged']:
        validate_run(medium_src)
        medium_out = medium_src
    else:
        medium_out = medium_dir / 'continuation-001'
        if not (medium_out / 'summary.json').exists():
            if medium_out.exists():
                raise RuntimeError(f'Incomplete continuation exists: {medium_out}')
            initial = medium_dir / 'initial-state-from-v12.json'
            if not initial.exists():
                state_file(medium_src, initial)
            config = read_json(medium_src / 'config.json')
            config['runtime']['max_iterations'] = MEDIUM_CONTINUATION
            config['output'].pop('directory', None)
            request = medium_dir / 'continuation-001.request.json'
            request.write_text(json.dumps(config, indent=2) + '\n')
            print(f'Continuing 256x64 for {MEDIUM_CONTINUATION} iterations from iteration {src_summary["iterations"]}', flush=True)
            code = run_cli(cli, request, medium_out, medium_dir / 'continuation-001.log', initial)
            (medium_out / 'provenance.json').write_text(json.dumps({
                'source_run': str(medium_src), 'source_iterations': src_summary['iterations'],
                'source_sha256': {k: digest(medium_src / k) for k in ('config.json', 'summary.json', 'final_state.vts')},
                'cli_sha256': digest(cli), 'source_commit': revision, 'exit_code': code}, indent=2) + '\n')
    s = read_json(medium_out / 'summary.json')
    medium_ok = s['converged']
    if medium_ok:
        validate_run(medium_out)
    manifest['runs'].append({'grid': [256, 64], 'directory': str(medium_out), 'converged': medium_ok, 'termination_reason': s['termination_reason'],
                             'source': 'v12' if medium_out == medium_src else 'v13 continuation', 'total_iterations': s['iterations'] + (0 if medium_out == medium_src else src_summary['iterations'])})
    (root / 'closure.json').write_text(json.dumps(manifest, indent=2) + '\n')
    if not medium_ok:
        print('Medium grid not closed under the V1.3 budget; fine grid not launched.', flush=True)
        return 0
    if args.skip_fine:
        return 0
    fine_dir = root / '512x128'
    fine_dir.mkdir(exist_ok=True)
    fine_out = fine_dir / 'attempt-001'
    if not (fine_out / 'summary.json').exists():
        if fine_out.exists():
            raise RuntimeError(f'Incomplete fine attempt exists: {fine_out}')
        config = read_json(coarse / 'config.json')
        config['mesh'] = {'nx': 512, 'nr': 128}
        config['runtime']['max_iterations'] = FINE_CAP
        config['output'].pop('directory', None)
        request = fine_dir / 'attempt-001.request.json'
        request.write_text(json.dumps(config, indent=2) + '\n')
        print(f'Running 512x128 with cap {FINE_CAP}', flush=True)
        code = run_cli(cli, request, fine_out, fine_dir / 'attempt-001.log')
        (fine_out / 'provenance.json').write_text(json.dumps({'cli_sha256': digest(cli), 'source_commit': revision, 'exit_code': code}, indent=2) + '\n')
    s = read_json(fine_out / 'summary.json')
    fine_ok = s['converged']
    if fine_ok:
        validate_run(fine_out)
    manifest['runs'].append({'grid': [512, 128], 'directory': str(fine_out), 'converged': fine_ok, 'termination_reason': s['termination_reason'], 'source': 'v13'})
    (root / 'closure.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (ValueError, RuntimeError, OSError) as error:
        print(f'V1.3 stopped: {error}', file=sys.stderr)
        sys.exit(2)

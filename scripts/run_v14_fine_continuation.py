#!/usr/bin/env python3
"""V1.4: continue the preserved V1.3 512x128 trajectory with the pre-registered budget.

See docs/verification/v14_fine_continuation_protocol.md. Writes runs/v14-closure/closure.json
listing the accepted V1.2 coarse run, the V1.3 medium continuation and the fine continuation.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from analyze_scientific_closure import read_json, validate_run
from run_v13_continuation import state_file, run_cli

ROOT = Path(__file__).resolve().parents[1]
FINE_CONTINUATION = 4500000


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--v13', type=Path, default=ROOT / 'runs/v13-closure')
    p.add_argument('--root', type=Path, default=ROOT / 'runs/v14-closure')
    p.add_argument('--cli', type=Path, default=ROOT / 'build/astraflow_cli')
    args = p.parse_args()
    v13, root, cli = args.v13.resolve(), args.root.resolve(), args.cli.resolve()
    root.mkdir(parents=True, exist_ok=True)
    manifest = read_json(v13 / 'closure.json')
    manifest['study'] = 'v14'
    accepted = [r for r in manifest['runs'] if r['converged']]
    for r in accepted:
        validate_run(Path(r['directory']))
    src = v13 / '512x128' / 'attempt-001'
    src_summary = read_json(src / 'summary.json')
    fine_dir = root / '512x128'
    fine_dir.mkdir(exist_ok=True)
    out = fine_dir / 'continuation-001'
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    if not (out / 'summary.json').exists():
        if out.exists():
            raise RuntimeError(f'Incomplete continuation exists: {out}')
        initial = fine_dir / 'initial-state-from-v13.json'
        if not initial.exists():
            state_file(src, initial)
        config = read_json(src / 'config.json')
        config['runtime']['max_iterations'] = FINE_CONTINUATION
        config['output'].pop('directory', None)
        request = fine_dir / 'continuation-001.request.json'
        request.write_text(json.dumps(config, indent=2) + '\n')
        print(f'Continuing 512x128 for {FINE_CONTINUATION} iterations from iteration {src_summary["iterations"]}', flush=True)
        code = run_cli(cli, request, out, fine_dir / 'continuation-001.log', initial)
        (out / 'provenance.json').write_text(json.dumps({
            'source_run': str(src), 'source_iterations': src_summary['iterations'],
            'source_sha256': {k: digest(src / k) for k in ('config.json', 'summary.json', 'final_state.vts')},
            'cli_sha256': digest(cli), 'source_commit': revision, 'exit_code': code}, indent=2) + '\n')
    s = read_json(out / 'summary.json')
    ok = s['converged']
    if ok:
        validate_run(out)
    manifest['runs'] = accepted + [{'grid': [512, 128], 'directory': str(out), 'converged': ok, 'termination_reason': s['termination_reason'],
                                    'source': 'v14 continuation', 'total_iterations': s['iterations'] + src_summary['iterations']}]
    (root / 'closure.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (ValueError, RuntimeError, OSError) as error:
        print(f'V1.4 stopped: {error}', file=sys.stderr)
        sys.exit(2)

#!/usr/bin/env python3
"""C2c guarded execution: coarse first, medium only after independent root-equivalence audit."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
from analyze_scientific_closure import read_json

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT/'runs/steady-acceleration-c2c'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_new(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'coarse', 'medium'))
    args = parser.parse_args()
    reference = ROOT/'runs/scientific-closure/128x32/attempt-001'
    if args.action == 'prepare':
        # Reproduction uses the same versioned bounds; never overwrite the recorded experiment.
        RUN.mkdir(parents=True, exist_ok=False)
        preset = read_json(ROOT/'docs/verification/steady_acceleration_c2c_preset.json')
        for name, expected in preset['coarse_source_sha256'].items():
            if digest(reference/name) != expected:
                raise ValueError('Accepted reference changed')
        write_new(RUN/'preset.json', preset)
        config = read_json(reference/'config.json')
        config['steady_solver'] = {'enabled':True, 'local_time_step':True, 'cfl':preset['steady_cfl']}
        config['output']['directory'] = ''
        write_new(RUN/'coarse-config.json', config)
        subprocess.run([str(ROOT/'build/astraflow_residual_probe'), 'balance', str(reference/'config.json'),
            str(ROOT/'runs/residual-closure-c2b/coarse-conservative.json'), str(RUN/'reference-balance.json')], check=True)
        return
    preset = read_json(RUN/'preset.json')
    for name, expected in preset['coarse_source_sha256'].items():
        if digest(reference/name) != expected:
            raise ValueError('Accepted reference changed')
    cli = ROOT/'build/astraflow_cli'
    if args.action == 'coarse':
        config = RUN/'coarse-config.json'
        out = RUN/'coarse-128x32'
        initial = None
    else:
        # Recompute from fields/history; do not trust a stored passing boolean.
        from analyze_steady_acceleration import coarse_evidence
        evidence = coarse_evidence(RUN)
        if not evidence['root_equivalence_pass']:
            raise RuntimeError('Coarse root-equivalence gate failed; medium is prohibited')
        out = RUN/'medium-256x64'
        if out.exists():
            raise FileExistsError('Preserving existing medium attempt')
        config = RUN/'medium-config.json'
        base = read_json(ROOT/'runs/residual-closure-c2b/warm-config.json')
        base['steady_solver'] = {'enabled': True, 'local_time_step': True, 'cfl': preset['steady_cfl']}
        base['runtime']['max_iterations'] = preset['medium_max_iterations']
        write_new(config, base)
        initial = ROOT/'runs/residual-closure-c2b/warm-initial.json'
    if out.exists():
        raise FileExistsError('Preserving existing attempt; no automatic reruns')
    cfg = read_json(config)
    original = read_json(reference/'config.json' if args.action == 'coarse' else
                         ROOT/'runs/residual-closure-c2b/warm-config.json')
    for key in original:
        if key != 'output' and original[key] != cfg[key]:
            raise ValueError(f'Production configuration changed: {key}')
    if cfg['steady_solver'] != {'enabled': True, 'local_time_step': True, 'cfl': preset['steady_cfl']}:
        raise ValueError('Steady settings differ from preset')
    execution = {'cli_sha256': digest(cli), 'configuration_sha256': digest(config),
                 'preset_sha256': digest(RUN/'preset.json'), 'time_integration': 'pseudo_time',
                 'initial_state_sha256': digest(initial) if initial else None,
                 'initial_condition': 'C2b conservative prolongation' if initial else 'original production initializer'}
    write_new(RUN/(args.action+'-execution.json'), execution)
    cmd = [str(cli), '--config', str(config), '--output', str(out)]
    if initial:
        cmd += ['--initial-state', str(initial)]
    subprocess.run(cmd, check=True)


if __name__ == '__main__':
    main()

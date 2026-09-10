#!/usr/bin/env python3
"""C2b preparation/execution; never reruns source closures or advances beyond 256x64."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
from analyze_scientific_closure import read_json, fields, prepare_diagnostic

ROOT = Path(__file__).resolve().parents[1]
SOURCE_FILES = ('config.json', 'summary.json', 'final_state.vts')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_new(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['prepare', 'probes', 'warm'])
    parser.add_argument('--root', type=Path, default=ROOT/'runs/residual-closure-c2b')
    parser.add_argument('--duration', type=float, default=20.)
    parser.add_argument('--sample-steps', type=int, default=1)
    args = parser.parse_args()
    root = args.root.resolve()
    root.mkdir(parents=True, exist_ok=True)
    coarse = ROOT/'runs/scientific-closure/128x32/attempt-001'
    medium = ROOT/'runs/scientific-closure/256x64/attempt-001'
    probe, cli = ROOT/'build/astraflow_residual_probe', ROOT/'build/astraflow_cli'

    if args.action == 'prepare':
        config = read_json(medium/'config.json')
        if read_json(coarse/'summary.json')['converged'] is not True:
            raise ValueError('Coarse source is not accepted')
        if [config['mesh']['nx'], config['mesh']['nr']] != [256, 64]:
            raise ValueError('Wrong saved medium grid')
        _, state = fields(coarse)
        write_new(root/'coarse-conservative.json', {
            'physical_conservative': list(map(list, zip(*state))),
            'provenance': {'source_directory': str(coarse),
                           'sha256': {k: digest(coarse/k) for k in SOURCE_FILES}}})
        config['output']['directory'] = ''
        write_new(root/'warm-config.json', config)
        subprocess.run([str(probe), 'prolong', str(coarse/'config.json'),
                        str(root/'coarse-conservative.json'), str(root/'warm-config.json'),
                        str(root/'warm-initial.json')], check=True)
        prepare_diagnostic(medium, root/'saved-medium')
        write_new(root/'sources.json', {
            'base_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT,
                                                    text=True).strip(),
            'medium_sha256': {k: digest(medium/k) for k in SOURCE_FILES},
            'coarse_sha256': {k: digest(coarse/k) for k in SOURCE_FILES},
            'probe_sha256': digest(probe), 'cli_sha256': digest(cli),
            'warm_initial_sha256': digest(root/'warm-initial.json')})
    elif args.action == 'probes':
        execution = {
            'probe_sha256': digest(probe),
            'normalized_config_sha256': digest(root/'saved-medium/normalized-config.json'),
            'primitives_sha256': digest(root/'saved-medium/primitives.txt'),
            'normalized_duration': args.duration, 'sampling_steps': args.sample_steps, 'cfl': .4}
        manifest = root/'probe-execution.json'
        if manifest.exists():
            if read_json(manifest) != execution:
                raise ValueError('Existing probe execution inputs/settings differ')
        else:
            write_new(manifest, execution)
        for mode in ('mc', 'van_leer', 'first_order'):
            out = root/f'probe-{mode}'
            if out.exists():
                summary = read_json(out/'summary.json')
                metadata = read_json(out/'metadata.json')
                if (summary['normalized_duration'] != args.duration or
                        metadata['sampling_steps'] != args.sample_steps or
                        metadata['mode'] != mode):
                    raise ValueError('Existing diagnostic settings mismatch')
                print(f'Preserving completed diagnostic {out}', flush=True)
                continue
            subprocess.run([str(probe), str(root/'saved-medium/normalized-config.json'),
                            str(root/'saved-medium/primitives.txt'), str(args.duration), '.4',
                            str(out), mode, str(args.sample_steps)], check=True)
    else:
        out = root/'warm-256x64'
        if out.exists():
            raise FileExistsError('Warm attempt exists; preserve it and choose a new explicit attempt')
        subprocess.run([str(cli), '--config', str(root/'warm-config.json'), '--initial-state',
                        str(root/'warm-initial.json'), '--output', str(out)], check=True)


if __name__ == '__main__':
    main()

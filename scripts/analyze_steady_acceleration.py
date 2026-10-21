#!/usr/bin/env python3
"""Independently audit the preset C2c root-equivalence gate and publish measured evidence."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
from analyze_scientific_closure import (read_json, read_csv, fields, norms, validate_run,
                                        OBSERVABLES, RESIDUALS)

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT/'runs/scientific-closure/128x32/attempt-001'
NAMES = ('mass', 'axial_momentum', 'radial_momentum', 'energy')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compact(summary):
    return {k:v for k,v in summary.items() if not k.startswith('station_')}


def compare_engineering(a, b):
    result = {}
    for key, value in a.items():
        if key not in b:
            continue
        other = b[key]
        if key.startswith('station_'):
            result[key] = {'physical':value, 'pseudo':other, 'norms':norms(value,other)}
        elif isinstance(value, (int,float)) and not isinstance(value,bool):
            result[key] = {'physical':value, 'pseudo':other,
                           'absolute_difference':abs(other-value),
                           'relative_difference':abs(other-value)/max(abs(value),1e-30)}
    return result


def coarse_evidence(root):
    root = Path(root)
    preset = read_json(root/'preset.json')
    for key, expected in preset['coarse_source_sha256'].items():
        if digest(REFERENCE/key) != expected:
            raise ValueError('Physical reference changed')
    path = root/'coarse-128x32'
    summary = read_json(path/'summary.json')
    reference = read_json(REFERENCE/'summary.json')
    config, reference_config = read_json(path/'config.json'), read_json(REFERENCE/'config.json')
    for key in reference_config:
        if key != 'output' and reference_config[key] != config[key]:
            raise ValueError(f'Production configuration changed: {key}')
    if config['steady_solver'] != {'enabled':True, 'local_time_step':True, 'cfl':preset['steady_cfl']}:
        raise ValueError('Unexpected steady settings')
    execution = read_json(root/'coarse-execution.json')
    if execution['preset_sha256'] != digest(root/'preset.json'):
        raise ValueError('Preset changed after experiment started')
    if execution['configuration_sha256'] != digest(root/'coarse-config.json'):
        raise ValueError('Input configuration changed after experiment started')
    result = {'preset':preset, 'execution':execution, 'source_hashes_verified':True,
              'physical':compact(reference), 'pseudo':compact(summary), 'gates':{},
              'root_equivalence_pass':False}
    gates = result['gates']
    try:
        validate_run(path)
        gates['unchanged_scientific_monitor'] = True
    except ValueError as error:
        gates['unchanged_scientific_monitor'] = False
        result['monitor_failure'] = str(error)
    history = read_csv(path/'convergence.csv')
    gates['physical_time_not_advanced'] = (summary['time_integration']=='pseudo_time' and
        summary['simulated_time']==0 and summary['physical_simulated_time']==0 and
        all(row['time']==0 and row['dt']==0 for row in history))
    gates['zero_fallbacks'] = summary['flux_fallbacks']==0 and summary['reconstruction_fallbacks']==0
    if not (path/'final_state.vts').exists():
        result['failure'] = 'No valid final field; do not advance to medium'
        return result
    profile, a = fields(REFERENCE)
    _, b = fields(path)
    scale = read_json(REFERENCE/'performance.json')['reference_scales']
    result['reference_scales'] = scale
    factors = [scale['density'],scale['density']*scale['velocity'],
               scale['density']*scale['velocity'],scale['pressure']]
    result['conservative_field_norms'] = {}
    result['field_maximum_locations'] = {}
    all_a, all_b = [], []
    for key, x, y, factor in zip(NAMES,a,b,factors):
        x,y = [v/factor for v in x], [v/factor for v in y]
        result['conservative_field_norms'][key] = norms(x,y)
        index = max(range(len(x)), key=lambda i:abs(x[i]-y[i]))
        nx = reference['grid'][0]
        result['field_maximum_locations'][key] = {'i':index%nx, 'j':index//nx,
            'x_over_L':profile['x_over_L'][index%nx], 'physical_normalized':x[index], 'pseudo_normalized':y[index]}
        all_a.extend(x); all_b.extend(y)
    result['conservative_field_norms']['all_components'] = norms(all_a,all_b)
    gates['field_linf'] = result['conservative_field_norms']['all_components']['Linf'] <= preset['field_linf']
    # Engineering fields come from the production analysis, excluding run/monitor metadata.
    engineering_keys = ('mass_flow','inlet_mass_flow','estimated_thrust','specific_impulse',
        'throat_mach','exit_mach','exit_pressure','exit_temperature','exit_axial_velocity',
        'exit_area','chamber_to_exit_pressure_ratio','mass_conservation_error','mass_flow_spread',
        'station_mass_flow','station_axial_mach')
    result['engineering'] = compare_engineering({k:reference[k] for k in engineering_keys},summary)
    gates['primary_engineering'] = all(result['engineering'][k]['relative_difference'] <=
        preset['primary_engineering_relative'] for k in OBSERVABLES)
    ref_balance, balance = read_json(root/'reference-balance.json'), read_json(path/'final_balance.json')
    result['physical_balance'], result['pseudo_balance'] = ref_balance, balance
    mismatch = balance['numerical_boundary_mass_flux_mismatch']
    gates['boundary_mass'] = (mismatch <= preset['boundary_mass_mismatch'] and
        abs(mismatch-ref_balance['numerical_boundary_mass_flux_mismatch']) <= preset['boundary_mass_mismatch_absolute_difference'])
    gates['completed_iterate_residual'] = max(balance['rms']) <= 1e-6
    result['physical_wall_seconds'] = read_json(REFERENCE/'performance.json')['wall_ms']/1000
    result['pseudo_wall_seconds'] = read_json(path/'performance.json')['wall_ms']/1000
    result['wall_speedup'] = result['physical_wall_seconds']/result['pseudo_wall_seconds']
    result['iteration_ratio_not_time_step_equivalence'] = reference['iterations']/summary['iterations']
    result['root_equivalence_pass'] = all(gates.values())
    result['failed_gates'] = [k for k,v in gates.items() if not v]
    return result


def plots(root, out, result):
    os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.cache/matplotlib'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    physical=read_csv(REFERENCE/'convergence.csv')
    pseudo=read_csv(root/'coarse-128x32/convergence.csv')
    fig,axes=plt.subplots(2,2,figsize=(10,6),layout='constrained')
    for ax,key in zip(axes.flat,RESIDUALS):
        for label,rows in [('Accepted physical SSP-RK2',physical),('Local pseudo RK2',pseudo)]:
            ax.semilogy([r['iteration'] for r in rows],[r[key] for r in rows],lw=1,label=label)
        ax.axhline(1e-6,color='black',ls='--',lw=.8)
        ax.set(xlabel='Iterations (different work and no equivalent elapsed time)',ylabel='Unscaled nondimensional RHS RMS',title=key.replace('_',' '))
        ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=8)
    fig.suptitle('128×32: physical second-stage RHS versus completed pseudo-iterate RHS\n'
                 + ('Preset field-equivalence gate PASSED' if result['root_equivalence_pass'] else
                    'Preset field-equivalence gate FAILED; medium prohibited'))
    fig.savefig(out/'c2c_coarse_residual_iterations.png',dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(6,4),layout='constrained')
    values=[result['physical_wall_seconds'],result['pseudo_wall_seconds']]
    bars=ax.bar(['Saved physical reference','Local pseudo-time'],values,color=['#2674a6','#dc8a24'])
    ax.bar_label(bars,labels=[f'{v:.3f} s' for v in values],padding=3)
    ax.set(ylabel='Total measured run wall time (s)',title=f'128×32 wall-time ratio: {result["wall_speedup"]:.2f}×\nTo monitor termination; separate equivalence gate '
           + ('passed' if result['root_equivalence_pass'] else 'FAILED'))
    ax.set_ylim(0,max(values)*1.15)
    fig.savefig(out/'c2c_coarse_wall_time.png',dpi=180);plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=ROOT/'runs/steady-acceleration-c2c')
    parser.add_argument('--output',type=Path,default=ROOT/'docs/verification')
    args=parser.parse_args()
    result=coarse_evidence(args.root)
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/'steady_acceleration_c2c_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    if 'wall_speedup' in result:
        plots(args.root,args.output,result)
    for name, path in (('physical_reference', REFERENCE), ('pseudo', args.root/'coarse-128x32')):
        (args.output/f'c2c_{name}_residual_history.csv').write_bytes((path/'convergence.csv').read_bytes())
    for key in ('gates','conservative_field_norms','physical_wall_seconds','pseudo_wall_seconds','wall_speedup','root_equivalence_pass'):
        print(key,result.get(key))


if __name__=='__main__':
    main()

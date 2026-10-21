#!/usr/bin/env python3
"""Analyze C2b saved-run evidence; never launches a solver or grants acceptance."""
import argparse
import csv
import gzip
import hashlib
import io
import os
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
COMPONENTS = ('mass', 'axial', 'radial', 'energy')
RESIDUALS = ('mass_residual', 'axial_momentum_residual',
             'radial_momentum_residual', 'energy_residual')
MODES = ('mc', 'van_leer', 'first_order')


def read_json(path):
    return json.loads(path.read_text())


def table(path):
    return np.genfromtxt(path, delimiter=',', names=True, encoding='utf-8')


def temporal(t, y):
    """Signed statistics; FFT only describes the sampled, linearly detrended signal."""
    slope, intercept = np.polyfit(t - t[0], y, 1)
    detrended = y - (intercept + slope * (t - t[0]))
    uniform = np.interp(np.linspace(t[0], t[-1], len(t)), t, detrended)
    power = np.abs(np.fft.rfft(uniform)) ** 2
    power[0] = 0
    if detrended.std() <= 32*np.finfo(float).eps*max(np.max(np.abs(y)), 1e-300):
        power[:] = 0
    peak = int(np.argmax(power))
    spacing = (t[-1] - t[0]) / (len(t) - 1)
    signs = np.sign(y[y != 0])
    return {'initial': float(y[0]), 'final': float(y[-1]),
            'minimum': float(y.min()), 'maximum': float(y.max()),
            'mean': float(y.mean()), 'stddev': float(y.std()),
            'linear_slope_per_normalized_time': float(slope),
            'detrended_stddev': float(detrended.std()),
            'sign_changes': int(np.count_nonzero(np.diff(signs))),
            'first_quarter_rms': float(np.sqrt(np.mean(y[:max(1, len(y)//4)]**2))),
            'last_quarter_rms': float(np.sqrt(np.mean(y[-max(1, len(y)//4):]**2))),
            'sampled_peak_period': float(len(t)*spacing/peak) if peak else None,
            'sampled_peak_power_fraction': float(power[peak]/power.sum()) if power.sum() else 0}


def cell_series(rows):
    branches = [n for n in rows.dtype.names if n.startswith('branch_')]
    result = {'i': int(rows['i'][0]), 'j': int(rows['j'][0]),
              'x_over_L': float(rows['x_over_L'][0]), 'samples': len(rows)}
    result['residual'] = {k: temporal(rows['time'], rows['rhs_'+k]) for k in COMPONENTS}
    result['branch_changes'] = {k: int(np.count_nonzero(np.diff(rows[k]))) for k in branches}
    result['branch_values'] = {k: [int(v) for v in np.unique(rows[k])] for k in branches}
    switched = np.zeros(len(rows)-1, dtype=bool)
    for key in branches:
        switched |= np.diff(rows[key]) != 0
    dy = np.abs(np.diff(rows['rhs_energy']))
    result['switching_intervals'] = int(switched.sum())
    result['mean_abs_energy_change_switch'] = float(dy[switched].mean()) if switched.any() else None
    result['mean_abs_energy_change_no_switch'] = float(dy[~switched].mean()) if (~switched).any() else None
    result['primitive_ranges'] = {k: {'initial': float(rows[k][0]), 'final': float(rows[k][-1]),
                                    'range': float(np.ptp(rows[k]))} for k in ('rho','u','v','p')}
    result['gradient_ranges'] = {k: [float(rows[k].min()), float(rows[k].max())]
                                 for k in rows.dtype.names if k.startswith(('gx_', 'gr_', 'sx_', 'sr_'))}
    result['maximum_reconstruction_faces'] = int(rows['reconstruction_faces'].max())
    result['maximum_fallback_faces'] = int(rows['fallback_faces'].max())
    return result


def snapshot(rows):
    result = {'top_cells': {}, 'component_rms': {}, 'first_four_rows_energy_square_fraction':
              float(np.sum(rows['rhs_energy'][rows['j'] < 4]**2)/np.sum(rows['rhs_energy']**2))}
    for key in COMPONENTS:
        result['component_rms'][key] = {p: float(np.sqrt(np.mean(rows[p+'_'+key]**2)))
                                        for p in ('rhs','convective','transport','source')}
        indices = np.argsort(np.abs(rows['rhs_'+key]))[-12:][::-1]
        names = ['cell','i','j','x_over_L','rho','u','v','p']
        names += [p+'_'+key for p in ('rhs','convective','transport','source')]
        result['top_cells'][key] = [{n: float(rows[n][i]) for n in names} for i in indices]
    return result


def compact_balance(balance):
    result = {k: v for k, v in balance.items() if k != 'interior_cell_estimate'}
    estimate = balance['interior_cell_estimate']
    result['interior_cell_estimate'] = {k: v for k, v in estimate.items() if not k.startswith('station_')}
    flux = balance['boundary_fluxes']
    inlet = -flux['inlet']['outward_total'][0]
    outlet = flux['outlet']['outward_total'][0]
    result['numerical_inlet_exit_mass_mismatch'] = abs(outlet-inlet)/max(abs(inlet),abs(outlet),1e-30)
    return result


def run_summary(path):
    s = read_json(path/'summary.json')
    return {k: v for k,v in s.items() if not k.startswith('station_')}


def convergence(path):
    rows = table(path/'convergence.csv')
    result = {'summary': run_summary(path), 'late_log_fits': {}}
    # Report finite-window fits, never extrapolated convergence as a measured result.
    for start in (50000, 150000, 200000):
        window = rows[rows['iteration'] >= start]
        if len(window) < 2:
            continue
        result['late_log_fits'][str(start)] = {}
        for key in RESIDUALS:
            fit = np.polyfit(window['time'], np.log(window[key]), 1)
            result['late_log_fits'][str(start)][key] = {
                'slope_per_second': float(fit[0]),
                'start': float(window[key][0]), 'end': float(window[key][-1]),
                'window_min': float(window[key].min()), 'window_max': float(window[key].max())}
    result['milestones'] = []
    for time in (.0001,.0005,.001,.002,.004,.006):
        if time <= rows['time'][-1]:
            result['milestones'].append({'time':time, **{k:float(np.interp(time,rows['time'],rows[k])) for k in RESIDUALS}})
    return result, rows


def analyze(root):
    result = {'schema':1, 'diagnostic_only_probes': True, 'sources':read_json(root/'sources.json'),
              'probe_execution':read_json(root/'probe-execution.json'),
              'prolongation':read_json(root/'warm-initial.json')['prolongation'], 'probes':{}}
    for kind, shape in (('coarse', '128x32'), ('medium', '256x64')):
        for name, expected in result['sources'][kind+'_sha256'].items():
            path = ROOT/'runs/scientific-closure'/shape/'attempt-001'/name
            assert hashlib.sha256(path.read_bytes()).hexdigest() == expected, path
    for name, key in (('normalized-config.json', 'normalized_config_sha256'),
                      ('primitives.txt', 'primitives_sha256')):
        assert hashlib.sha256((root/'saved-medium'/name).read_bytes()).hexdigest() == result['probe_execution'][key]
    assert hashlib.sha256((root/'warm-initial.json').read_bytes()).hexdigest() == result['sources']['warm_initial_sha256']
    result['source_hashes_verified'] = True
    scale = read_json(ROOT/'runs/scientific-closure/256x64/attempt-001/performance.json')['reference_scales']
    result['reference_scales'] = scale
    raw = {}
    common = None
    for mode in MODES:
        path = root/f'probe-{mode}'
        metadata = read_json(path/'metadata.json')
        if common is None:
            common = metadata['watch_cells']
        assert common == metadata['watch_cells']
        balances = [json.loads(s) for s in (path/'balance.jsonl').read_text().splitlines()]
        tracks = table(path/'tracked.csv')
        initial, final = table(path/'initial.csv'), table(path/'final.csv')
        assert np.all(np.isfinite(tracks.view(float)))
        assert balances[-1]['time'] == metadata['normalized_duration']
        data = {'metadata':metadata, 'physical_duration_seconds':balances[-1]['time']*scale['time'],
                'initial_balance':compact_balance(balances[0]), 'final_balance':compact_balance(balances[-1]),
                'initial_snapshot':snapshot(initial), 'final_snapshot':snapshot(final),
                'max_relative_assembly_error':np.max([b['relative_identity_error'] for b in balances],axis=0).tolist(),
                'tracked_cells':{str(int(i)):cell_series(tracks[tracks['cell']==i]) for i in np.unique(tracks['cell'])}}
        data['global_rms_history'] = {k:temporal(np.array([b['time'] for b in balances]),
            np.array([b['rms'][j] for b in balances])) for j,k in enumerate(COMPONENTS)}
        result['probes'][mode] = data
        raw[mode] = (balances, tracks, initial, final)
    # Ensure all three experiments actually start with identical primitive data.
    for mode in MODES[1:]:
        for key in ('rho','u','v','p'):
            assert np.array_equal(raw['mc'][2][key],raw[mode][2][key])
    result['identical_probe_initial_primitives'] = True
    result['original_medium'], cold = convergence(ROOT/'runs/scientific-closure/256x64/attempt-001')
    result['warm_medium'], warm = convergence(root/'warm-256x64')
    baseline = read_json(ROOT/'runs/scientific-closure/256x64/attempt-001/config.json')
    warm_config = read_json(root/'warm-256x64/config.json')
    for key in baseline:
        if key != 'output':
            assert baseline[key] == warm_config[key], key
    result['unchanged_warm_production_configuration'] = True
    return result, raw, cold, warm


def plots(out, raw, cold, warm, scale):
    os.environ.setdefault('MPLCONFIGDIR', str(ROOT/'.cache/matplotlib'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':9, 'axes.grid':True, 'grid.alpha':.2})
    fig, axes = plt.subplots(2,2,figsize=(10,6),layout='constrained')
    for ax, key in zip(axes.flat, RESIDUALS):
        for label, rows in [('Original saved cold run',cold),('Conservative coarse initialization',warm)]:
            ax.semilogy(rows['time']*1e3,rows[key],label=label,lw=.8)
        ax.axhline(1e-6,color='black',ls='--',lw=.8)
        ax.set(xlabel='Physical time since initialization (ms)',ylabel='Nondimensional stage RHS RMS',title=key.replace('_',' '))
    axes[0,0].legend(fontsize=7)
    fig.savefig(out/'c2b_warm_comparison.png',dpi=180); plt.close(fig)
    fig, axes = plt.subplots(2,2,figsize=(10,6),layout='constrained')
    for ax, key, k in zip(axes.flat,COMPONENTS,range(4)):
        for mode, (balance,*_) in raw.items():
            ax.semilogy([b['time']*scale*1e6 for b in balance],[b['rms'][k] for b in balance],label=mode)
        ax.set(xlabel='Continuation physical time (µs)',ylabel='Snapshot RHS RMS',title=key)
    axes[0,0].legend()
    fig.suptitle('Equal-time diagnostics from the same saved medium state; altered schemes are ineligible')
    fig.savefig(out/'c2b_equal_time.png',dpi=180); plt.close(fig)
    fig, axes = plt.subplots(3,2,figsize=(10,8),layout='constrained')
    for col, cell in enumerate((38,110)):
        tracks=raw['mc'][1]; rows=tracks[(tracks['cell']==cell) & (tracks['time'] <= 2)]; t=rows['time']*scale*1e6
        axes[0,col].plot(t,rows['rhs_energy'],lw=.7)
        axes[0,col].set(title=f'MC cell (i={cell}, j=0), x/L={rows["x_over_L"][0]:.4f}',ylabel='Signed energy RHS')
        for key in ('branch_x_rho','branch_x_u','branch_r_rho','branch_r_v'):
            axes[1,col].step(t,rows[key],where='post',lw=.7,label=key)
        axes[1,col].set(ylabel='Limiter branch code'); axes[1,col].legend(fontsize=6)
        for key in ('rho','u','p'):
            axes[2,col].plot(t,(rows[key]-rows[key][0])/abs(rows[key][0]),lw=.7,label=key)
        axes[2,col].set(xlabel='Continuation physical time (µs)',ylabel='Relative primitive change')
        axes[2,col].legend(fontsize=7)
    fig.suptitle('MC signed oscillations: first 10 µs of the 100.4 µs continuation')
    fig.savefig(out/'c2b_tracked_cells.png',dpi=180); plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=ROOT/'runs/residual-closure-c2b')
    p.add_argument('--output',type=Path,default=ROOT/'docs/verification')
    args=p.parse_args(); args.output.mkdir(parents=True,exist_ok=True)
    result, raw, cold, warm = analyze(args.root)
    (args.output/'residual_closure_c2b_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    plots(args.output,raw,cold,warm,result['reference_scales']['time'])
    # Publish lossless, full-cadence focus traces and global RMS histories.
    stream = io.StringIO(newline='')
    writer = csv.writer(stream, lineterminator='\n')
    names = raw['mc'][1].dtype.names
    writer.writerow(('mode', *names))
    for mode, (_, tracks, *_) in raw.items():
        for row in tracks[np.isin(tracks['cell'], (38,110))]:
            writer.writerow((mode, *(format(float(row[n]), '.17g') for n in names)))
    (args.output/'c2b_focus_traces.csv.gz').write_bytes(gzip.compress(stream.getvalue().encode(), mtime=0))
    with (args.output/'c2b_rms_history.csv').open('w', newline='') as stream:
        writer = csv.writer(stream, lineterminator='\n'); writer.writerow(('mode','iteration','normalized_time',*COMPONENTS))
        for mode,(balance,*_) in raw.items():
            for b in balance:
                writer.writerow((mode,b['iteration'],b['time'],*b['rms']))
    for mode, data in result['probes'].items():
        print(mode, data['initial_balance']['rms'], '->',data['final_balance']['rms'])
    print('warm',result['warm_medium']['summary']['termination'],result['warm_medium']['summary']['final_residuals'])


if __name__ == '__main__':
    main()

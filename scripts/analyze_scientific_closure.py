#!/usr/bin/env python3
"""Independent closure statistics and three-grid Richardson/GCI analysis (stdlib only)."""
import argparse
import csv
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET

OBSERVABLES = ('mass_flow', 'estimated_thrust', 'specific_impulse', 'throat_mach', 'exit_mach', 'exit_pressure')
GOALS = dict(zip(OBSERVABLES, (0.005, 0.005, 0.005, 0.005, 0.005, 0.01)))
METRICS = OBSERVABLES + ('exit_axial_velocity', 'exit_temperature', 'mass_conservation_error', 'mass_flow_spread')
RESIDUALS = ('mass_residual', 'axial_momentum_residual', 'radial_momentum_residual', 'energy_residual')

def read_json(path):
    def invalid(value):
        raise ValueError(f'Nonfinite JSON constant: {value}')
    return json.loads(Path(path).read_text(), parse_constant=invalid)

def read_csv(path):
    with Path(path).open() as source:
        rows = [{k: float(v) for k, v in row.items()} for row in csv.DictReader(source)]
    if not rows or any(not math.isfinite(x) for row in rows for x in row.values()):
        raise ValueError(f'Empty or nonfinite history: {path}')
    if any(b['iteration'] <= a['iteration'] for a, b in zip(rows, rows[1:])):
        raise ValueError(f'Unordered history: {path}')
    return rows

def stability(x, y):
    if len(x) != len(y) or len(x) < 2 or any(not math.isfinite(v) for v in [*x, *y]):
        raise ValueError('Stability requires at least two finite samples')
    mean = math.fsum(y) / len(y)
    xm = math.fsum(x) / len(x)
    scale = max(abs(mean), 1e-30)
    xx = math.fsum((v-xm)**2 for v in x)
    if xx == 0:
        raise ValueError('Zero observation interval')
    slope = math.fsum((a-xm)*(b-mean) for a,b in zip(x,y)) / xx
    return {'mean': mean, 'relative_range': (max(y)-min(y))/scale,
            'relative_stddev': math.sqrt(math.fsum((v-mean)**2 for v in y)/len(y))/scale,
            'relative_drift': slope*(max(x)-min(x))/scale}

def richardson(coarse, medium, fine, r=2.0, safety=1.25):
    result = {'valid': False, 'order': None, 'extrapolated': None, 'gci_percent': None}
    if not all(math.isfinite(x) for x in (coarse,medium,fine,r,safety)) or r <= 1 or safety < 1:
        raise ValueError('Invalid Richardson inputs')
    a, b = coarse-medium, medium-fine
    floor = 32 * math.ulp(max(abs(coarse),abs(medium),abs(fine)))
    if min(abs(a),abs(b)) <= floor:
        result['reason'] = 'differences indistinguishable from roundoff'
    elif a*b <= 0:
        result['reason'] = 'oscillatory / non-monotonic sequence'
    elif abs(a) <= abs(b):
        result['reason'] = 'differences do not decrease under refinement'
    elif fine == 0:
        result['reason'] = 'relative GCI undefined at zero fine value'
    else:
        p = math.log(abs(a/b)) / math.log(r)
        denominator = math.expm1(p*math.log(r))
        result.update(valid=True, order=p, extrapolated=fine+(fine-medium)/denominator,
                      gci_percent=100*safety*abs((fine-medium)/fine)/denominator,
                      reason='monotonic decreasing differences; conditional single-power asymptotic model',
                      safety_factor=safety, refinement_ratio=r)
    return result

def norms(a,b):
    if len(a)!=len(b) or not a:
        raise ValueError('Norm arrays must have equal positive length')
    d=[abs(x-y) for x,y in zip(a,b)]
    if not all(math.isfinite(x) for x in d): raise ValueError('Nonfinite norm data')
    return {'L1':math.fsum(d)/len(d), 'L2':math.sqrt(math.fsum(x*x for x in d)/len(d)), 'Linf':max(d)}

def validate_run(directory):
    directory=Path(directory)
    s=read_json(directory/'summary.json'); config=read_json(directory/'config.json')
    history=read_csv(directory/'engineering.csv'); residual=read_csv(directory/'convergence.csv')
    if s.get('termination_reason')!='steady_converged' or s.get('converged') is not True:
        raise ValueError(f'{directory}: not steady-converged ({s.get("termination_reason")})')
    if s['grid'] != [config['mesh']['nx'],config['mesh']['nr']] or any(s[k]!=config['runtime'][k] for k in ('backend','precision')):
        raise ValueError('Configuration/summary identity mismatch')
    end=s['iterations']; settings=config['convergence']; start=end-settings['window_iterations']
    if end < settings['minimum_iterations'] or history[-1]['iteration'] != end or residual[-1]['iteration'] != end:
        raise ValueError('Incomplete final history or minimum iteration count')
    window=[row for row in history if row['iteration']>=start]
    if window[0]['iteration']!=start or window[-1]['iteration']!=end:
        raise ValueError('Incomplete engineering window')
    if any(b['iteration']-a['iteration']!=settings['sampling_interval'] for a,b in zip(window,window[1:])):
        raise ValueError('Missing engineering samples')
    rw=[row for row in residual if row['iteration']>=start]
    if not {row['iteration'] for row in window}.issubset({row['iteration'] for row in rw}):
        raise ValueError('Missing sampled residuals')
    target=config['runtime']['residual_tolerance']
    if not 0 < target <= 1e-6 or any(row[k]>target for row in rw for k in RESIDUALS):
        raise ValueError('Absolute residual closure failed')
    relative=settings['relative_residual_target']
    if relative and any(row[k]>s['initial_residuals'][i]*relative for row in rw for i,k in enumerate(RESIDUALS)):
        raise ValueError('Relative residual closure failed')
    st={k:stability([r['iteration'] for r in window],[r[k] for r in window]) for k in OBSERVABLES}
    if any(v['relative_range']>min(5e-4,settings['observable_tolerance']) for v in st.values()):
        raise ValueError('Engineering stability failed')
    for k,limit,setting in [('mass_conservation_error',1e-3,'mass_mismatch_tolerance'),('mass_flow_spread',2e-3,'station_spread_tolerance')]:
        if any(row[k]>min(limit,settings[setting]) for row in window): raise ValueError(f'{k} failed')
    for i,k in enumerate(RESIDUALS):
        if not math.isclose(residual[-1][k],s['final_residuals'][i],rel_tol=1e-12,abs_tol=1e-30):
            raise ValueError('Summary residual mismatch')
    for k in METRICS:
        if not math.isfinite(s[k]): raise ValueError(f'Nonfinite engineering metric {k}')
    for k in OBSERVABLES + ('mass_conservation_error','mass_flow_spread'):
        if not math.isclose(history[-1][k],s[k],rel_tol=1e-12,abs_tol=1e-30):
            raise ValueError('Summary engineering mismatch')
    if not (directory/'final_state.vts').exists(): raise ValueError('Missing final field')
    return s,st

def fields(directory):
    directory=Path(directory); s=read_json(directory/'summary.json'); nx,nr=s['grid']
    root=ET.parse(directory/'final_state.vts')
    piece=root.find('./StructuredGrid/Piece')
    arrays={a.attrib['Name']:[float(v) for v in a.text.split()] for a in piece.findall('./CellData/DataArray')}
    for k,v in arrays.items():
        if len(v)!=nx*nr*(3 if k=='velocity' else 1) or not all(math.isfinite(x) for x in v): raise ValueError('Invalid VTK array')
    for k in ('density','pressure','temperature'):
        if min(arrays[k])<=0: raise ValueError('Nonpositive VTK state')
    points=[float(x) for x in piece.find('./Points/DataArray').text.split()]
    x=[(points[3*i]+points[3*(i+1)])/2 for i in range(nx)]
    length=points[3*nx]-points[0]
    profiles={'x_over_L':[(v-points[0])/length for v in x], 'centreline_mach':arrays['Mach'][:nx],
              'centreline_pressure':arrays['pressure'][:nx], 'centreline_temperature':arrays['temperature'][:nx],
              'wall_pressure':arrays['pressure'][(nr-1)*nx:]}
    rho=arrays['density']; vel=arrays['velocity']; energy=arrays['total_energy']
    conservative=[rho,[rho[i]*vel[3*i] for i in range(nx*nr)], [rho[i]*vel[3*i+1] for i in range(nx*nr)], [rho[i]*energy[i] for i in range(nx*nr)]]
    return profiles,conservative

def interpolate(x,y,target):
    import bisect
    result=[]
    for t in target:
        if not x[0]<=t<=x[-1]: raise ValueError('Profile extrapolation prohibited')
        i=max(0,min(len(x)-2,bisect.bisect_right(x,t)-1))
        result.append(y[i]+(y[i+1]-y[i])*(t-x[i])/(x[i+1]-x[i]))
    return result

def analyze(root):
    root=Path(root); manifest=read_json(root/'closure.json')
    runs=[(root/e['directory']).resolve() for e in manifest['runs'] if e.get('converged')]
    verified=[validate_run(p)[0] for p in runs]
    if len(verified)<1: raise ValueError('No accepted steady grids')
    def physical_signature(p):
        config=read_json(p/'config.json')
        config.pop('mesh');config.pop('output')
        config['runtime'].pop('max_iterations')
        return config
    if any(physical_signature(p)!=physical_signature(runs[0]) for p in runs[1:]):
        raise ValueError('Grid study mixes physical/numerical configuration or precision')
    profiles=[fields(p)[0] for p in runs]
    (root/'profiles').mkdir(exist_ok=True)
    for s,p in zip(verified,profiles):
        with (root/'profiles'/f'{s["grid"][0]}x{s["grid"][1]}.csv').open('w') as f:
            w=csv.writer(f);w.writerow(p);w.writerows(zip(*p.values()))
    comparisons=[]
    for coarse,fine,a,b in zip(verified,verified[1:],profiles,profiles[1:]):
        if any(f != 2*c for c,f in zip(coarse['grid'],fine['grid'])): raise ValueError('Expected r=2 in both directions')
        changes={k:abs((fine[k]-coarse[k])/fine[k]) for k in METRICS if fine[k]!=0}
        lo=max(a['x_over_L'][0],b['x_over_L'][0]); hi=min(a['x_over_L'][-1],b['x_over_L'][-1])
        x=[lo+(hi-lo)*i/2048 for i in range(2049)]
        pn={k:norms(interpolate(a['x_over_L'],a[k],x),interpolate(b['x_over_L'],b[k],x)) for k in a if k!='x_over_L'}
        geometry=read_json(runs[0]/'config.json')['geometry']
        throat=(geometry['chamber_length']+geometry['contraction_length'])/(geometry['chamber_length']+geometry['contraction_length']+geometry['expansion_length'])
        throat_indices=[i for i,t in enumerate(x) if abs(t-throat)<=.03]
        for k in pn:
            aa=interpolate(a['x_over_L'],a[k],x);bb=interpolate(b['x_over_L'],b[k],x)
            pn[k]['Linf_x_over_L']=x[max(range(len(x)),key=lambda i:abs(aa[i]-bb[i]))]
            pn[k]['throat_band_norms']=norms([aa[i] for i in throat_indices],[bb[i] for i in throat_indices])
        comparisons.append({'coarse':coarse['grid'],'fine':fine['grid'],'relative_differences':changes,'profile_norms':pn,
            'common_x_over_L':[lo,hi],'common_points':len(x),'global_goals_pass':all(changes[k]<=v for k,v in GOALS.items())})
    rich={}
    if len(verified)>=3:
        for k in OBSERVABLES:
            rich[k]=richardson(*(s[k] for s in verified[-3:]))
    report={'excluded_runs':[e for e in manifest['runs'] if not e.get('converged')], 'grids':[{k:s[k] for k in ('grid','iterations','simulated_time',*METRICS)} for s in verified],
            'comparisons':comparisons,'richardson':rich,
            'global_grid_goals_pass':len(verified)>=3 and comparisons[-1]['global_goals_pass'],
            'valid_gci_goals_pass':bool(rich) and all(not v['valid'] or v['gci_percent']<=1 for v in rich.values()),
            'gci_assumption':'Conditional on a single leading power of h. Three samples alone do not prove the asymptotic regime.',
            'profile_definition':'First radial cell row approximates centreline; last radial row approximates wall pressure. No extrapolation to boundaries; common axial overlap includes the throat.'}
    (root/'analysis.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    with (root/'grid_metrics.csv').open('w') as f:
        w=csv.writer(f); w.writerow(('nx','nr','iterations',*METRICS))
        for s in verified: w.writerow((*s['grid'],s['iterations'],*(s[k] for k in METRICS)))
    return report

def compare_backends(cpu_directory, cuda_directory):
    cpu_directory, cuda_directory = Path(cpu_directory), Path(cuda_directory)
    a,_=validate_run(cpu_directory);b,_=validate_run(cuda_directory)
    if a['grid']!=b['grid'] or a['precision']!='double' or b['precision']!='double':
        raise ValueError('CPU/CUDA reference comparison requires matching FP64 grids')
    def signature(directory):
        c=read_json(directory/'config.json');c.pop('output')
        c['runtime'].pop('backend');c['runtime'].pop('max_iterations')
        return c
    if signature(cpu_directory)!=signature(cuda_directory):
        raise ValueError('CPU/CUDA comparison mixes configuration')
    _,u=fields(cpu_directory);_,v=fields(cuda_directory)
    ref=read_json(cpu_directory/'performance.json')['reference_scales']
    factors=[ref['density'],ref['density']*ref['velocity'],ref['density']*ref['velocity'],ref['pressure']]
    field_norms={name:norms([x/f for x in aa],[x/f for x in bb]) for name,aa,bb,f in zip(('rho','rho_u','rho_v','rho_E'),u,v,factors)}
    engineering={k:abs((a[k]-b[k])/b[k]) for k in METRICS if b[k]!=0}
    return {'cpu_iterations':a['iterations'],'cuda_iterations':b['iterations'],
            'nondimensional_conservative_norms':field_norms,'engineering_relative_differences':engineering,
            'maximum_field_difference':max(x['Linf'] for x in field_norms.values()),
            'same_iteration_parity_bound':1e-10,
            'parity_pass':a['iterations']==b['iterations'] and max(x['Linf'] for x in field_norms.values())<=1e-10}

def prepare_diagnostic(directory, output):
    directory,output=Path(directory),Path(output)
    output.mkdir(parents=True,exist_ok=True)
    config_path,primitive_path=output/'normalized-config.json',output/'primitives.txt'
    if config_path.exists() or primitive_path.exists():
        raise FileExistsError('Diagnostic inputs already exist; choose a new output directory')
    c=read_json(directory/'config.json')
    scale=read_json(directory/'performance.json')['reference_scales']
    length,rho,velocity,pressure,temperature=(scale[k] for k in ('length','density','velocity','pressure','temperature'))
    for key in c['geometry']: c['geometry'][key]/=length
    c['gas']['gas_constant']=1
    c['gas']['viscosity']/=rho*velocity*length
    c['gas']['rho_floor']/=rho;c['gas']['p_floor']/=pressure;c['gas']['temperature_floor']/=temperature
    for key,factor in [('P0',pressure),('T0',temperature),('back_pressure',pressure),('wall_temperature',temperature),('upper_wall_speed',velocity)]:
        c['boundary_conditions'][key]/=factor
    _,u=fields(directory)
    config_path.write_text(json.dumps(c,indent=2,allow_nan=False)+'\n')
    with primitive_path.open('x') as f:
        for d,mx,mr,e in zip(*u):
            p=(c['gas']['gamma']-1)*(e-.5*(mx*mx+mr*mr)/d)
            f.write(f'{d/rho:.17g} {mx/d/velocity:.17g} {mr/d/velocity:.17g} {p/pressure:.17g}\n')
    return {'configuration':str(config_path),'primitives':str(primitive_path)}

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('root',type=Path)
    parser.add_argument('--prepare-diagnostic',type=Path,metavar='NEW_DIRECTORY')
    parser.add_argument('--compare-cpu',type=Path,metavar='CPU_RUN')
    args=parser.parse_args()
    if args.prepare_diagnostic:
        result=prepare_diagnostic(args.root,args.prepare_diagnostic)
    elif args.compare_cpu:
        result=compare_backends(args.compare_cpu,args.root)
    else:
        result=analyze(args.root)
    print(json.dumps(result,indent=2,allow_nan=False))

#!/usr/bin/env python3
"""Sequential, resumable CLI-driven closure study; completed outputs are immutable."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time
from analyze_scientific_closure import read_json, validate_run

ROOT=Path(__file__).resolve().parents[1]
def digest(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def atomic_json(path,value):
    tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');tmp.replace(path)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=ROOT/'runs/scientific-closure')
    p.add_argument('--config',type=Path,default=ROOT/'examples/rocket_nozzle/scientific_closure.json')
    p.add_argument('--cli',type=Path,default=ROOT/'build/astraflow_cli')
    p.add_argument('--backend',choices=['cpu','cuda'],default='cuda')
    p.add_argument('--precision',choices=['float','double'],default='double')
    p.add_argument('--iteration-limit', action='append', default=[], metavar='GRID=COUNT')
    p.add_argument('--grids',nargs='+',default=['128x32','256x64','512x128'])
    args=p.parse_args();root=args.root.resolve();root.mkdir(parents=True,exist_ok=True)
    limits={key:int(value) for key,value in (item.split('=') for item in args.iteration_limit)}
    base=read_json(args.config);cli=args.cli.resolve();records=[]
    manifest={'schema':1,'backend':args.backend,'precision':args.precision,'runs':records}
    revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    # Exclusive lock prevents concurrent orchestrators from claiming the same attempt directory.
    import fcntl
    with (root/'.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        for grid in args.grids:
            nx,nr=map(int,grid.split('x'))
            config=json.loads(json.dumps(base));config['mesh']={'nx':nx,'nr':nr}
            config['runtime'].update(backend=args.backend,precision=args.precision)
            if grid in limits: config['runtime']['max_iterations']=limits[grid]
            config['output'].pop('directory',None)
            provenance={'config':config,'cli_sha256':digest(cli)}
            parent=root/grid;parent.mkdir(exist_ok=True)
            previous=sorted(parent.glob('attempt-*'))
            reused=None
            for attempt in previous:
                try:
                    meta=read_json(attempt/'provenance.json')
                    if meta['request']!=provenance: continue
                    if meta['exit_code']!=0: raise ValueError('Previous CLI process failed')
                    if any(digest(attempt/k)!=v for k,v in meta['sha256'].items()): raise ValueError('Output checksum mismatch')
                    summary,_=validate_run(attempt)
                    reused=(attempt,summary);break
                except (ValueError,KeyError,OSError) as error:
                    print(f'Preserving incomplete/incompatible attempt {attempt}: {error}',flush=True)
            if reused:
                attempt,summary=reused
                print(f'Reusing independently verified {attempt}',flush=True)
            else:
                attempt=parent/f'attempt-{len(previous)+1:03d}'
                request=parent/f'request-{len(previous)+1:03d}.json'
                if request.exists(): raise RuntimeError(f'Request path already exists: {request}')
                attempt.mkdir()  # Preserve a failed startup as a distinct attempt too.
                request.write_text(json.dumps(config,indent=2)+'\n')
                start=time.time()
                print(f'Running {grid} {args.backend} {args.precision}: {attempt}',flush=True)
                with (parent/f'cli-{len(previous)+1:03d}.log').open('x') as log:
                    process=subprocess.run([str(cli),'--config',str(request),'--output',str(attempt)],stdout=log,stderr=subprocess.STDOUT)
                summary=read_json(attempt/'summary.json') if (attempt/'summary.json').exists() else {'converged':False,'termination_reason':'incomplete'}
                meta={'request':provenance,'source_commit':revision,'exit_code':process.returncode,'wall_seconds':time.time()-start,
                      'sha256':{k:digest(attempt/k) for k in ('config.json','summary.json','performance.json','engineering.csv','convergence.csv','final_state.vts') if (attempt/k).exists()}}
                atomic_json(attempt/'provenance.json',meta)
            records.append({'grid':[nx,nr],'directory':str(attempt.relative_to(root)),
                            'converged':summary.get('converged',False),'termination_reason':summary.get('termination_reason')})
            atomic_json(root/'closure.json',manifest)
            validate_run(attempt)  # Never proceed to a finer grid after failed closure.
    return 0
if __name__=='__main__':
    try: sys.exit(main())
    except (ValueError,RuntimeError,OSError) as error:
        print(f'Closure stopped: {error}',file=sys.stderr);sys.exit(2)

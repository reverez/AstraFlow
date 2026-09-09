#!/usr/bin/env python3
"""Generate figures directly from preserved CFD outputs; Agg, no interactive display."""
import argparse
from pathlib import Path
import os
os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[1] / ".cache/matplotlib"))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyze_scientific_closure import read_json, read_csv, fields, RESIDUALS

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root',type=Path)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--history-grid', help='Select a diagnostic history, e.g. 256x64; termination remains labelled')
    args=parser.parse_args();root=args.root;output=args.output or root/'figures';output.mkdir(parents=True,exist_ok=True)
    manifest=read_json(root/'closure.json')
    runs=[root/e['directory'] for e in manifest['runs']]
    summaries=[read_json(p/'summary.json') for p in runs]
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                         'savefig.dpi':180,'figure.constrained_layout.use':True})
    def save(fig,name):
        fig.savefig(output/name,facecolor='white');plt.close(fig)
    # Always identify termination; an iteration-limited diagnostic plot cannot imply closure.
    selected=[(p,s) for p,s in zip(runs,summaries) if s['converged']]
    if args.history_grid:
        selected=[(p,s) for p,s in zip(runs,summaries) if f'{s["grid"][0]}x{s["grid"][1]}'==args.history_grid]
        if not selected: raise ValueError('Requested history grid not in manifest')
    p,s=(selected or list(zip(runs,summaries)))[-1];label=f'{s["grid"][0]} × {s["grid"][1]}, {s["backend"]} {s["precision"]}, {s["termination_reason"]}'
    rows=read_csv(p/'convergence.csv');x=[r['iteration'] for r in rows]
    fig,ax=plt.subplots(figsize=(8,4.4))
    for k,name in zip(RESIDUALS,('Mass','Axial momentum','Radial momentum','Energy')):
        ax.semilogy(x,[r[k] for r in rows],label=name,lw=1.1)
    ax.axhline(1e-6,color='black',ls='--',lw=1,label='Absolute target')
    ax.set(xlabel='Iteration',ylabel='Nondimensional RMS time derivative',title=label,xlim=(0,x[-1]))
    ax.grid(alpha=.2);ax.legend(ncol=2,fontsize=9);save(fig,'rocket-residual-convergence.png')
    rows=read_csv(p/'engineering.csv');x=[r['iteration'] for r in rows]
    fig,axes=plt.subplots(3,1,figsize=(8,7),sharex=True)
    for ax,k,name in zip(axes,('mass_flow','estimated_thrust','exit_mach'),('Mass flow [kg/s]','Thrust [N]','Exit Mach')):
        ax.plot(x,[r[k] for r in rows],lw=1.2);ax.set_ylabel(name);ax.grid(alpha=.2)
        ax.ticklabel_format(axis='y',style='plain',useOffset=False)
    axes[0].set_title(label+'\nFull recorded range; axes resolve transient changes')
    axes[-1].set(xlabel='Iteration',xlim=(0,x[-1]));save(fig,'rocket-engineering-convergence.png')
    accepted=[(p,s) for p,s in zip(runs,summaries) if s['converged']]
    if len(accepted)>=2:
        analysis=read_json(root/'analysis.json')
        fig,axes=plt.subplots(2,3,figsize=(11,6),sharex=True)
        for ax,k,name in zip(axes.flat,('mass_flow','estimated_thrust','specific_impulse','throat_mach','exit_mach','exit_pressure'),('Mass flow [kg/s]','Thrust [N]','Isp [s]','Throat Mach','Exit Mach','Exit pressure [kPa]')):
            factor=.001 if k=='exit_pressure' else 1
            ax.plot([1/s['grid'][0] for _,s in accepted],[factor*s[k] for _,s in accepted],'-o',lw=1.1,ms=4)
            ax.set_ylabel(name);ax.grid(alpha=.2);ax.ticklabel_format(axis='y',style='plain',useOffset=False)
        for ax in axes[-1]:ax.set_xlabel('Axial spacing measure, 1 / nx')
        fig.suptitle('Steady rocket grid study — systematic refinement\nAxes resolve measured discretisation changes')
        save(fig,'rocket-grid-convergence.png')
    for name,key,ylabel,scale in [('rocket-centreline-mach.png','centreline_mach','Centreline-row Mach',1),('rocket-wall-pressure.png','wall_pressure','Wall-row pressure [kPa]',.001)]:
        fig,ax=plt.subplots(figsize=(8,4.4))
        for p,s in accepted:
            profile,_=fields(p);ax.plot(profile['x_over_L'],[scale*v for v in profile[key]],label=f'{s["grid"][0]} × {s["grid"][1]}',lw=1.2)
        ax.set(xlabel='Axial position, x/L',ylabel=ylabel,xlim=(0,1));ax.grid(alpha=.2);ax.legend()
        save(fig,name)
    print(f'Wrote actual-data figures to {output}; accepted grids: {len(accepted)}. Grid plot requires at least two converged grids.')
if __name__=='__main__':main()

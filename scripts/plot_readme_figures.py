#!/usr/bin/env python3
"""README figures rendered from preserved AstraFlow outputs only (Agg backend, no display).

Subcommands:
  sod RUN_DIR            numerical Sod solution against an independent exact Riemann solution
  field RUN_DIR          Mach and pressure fields of a final_state.vts, mirrored about the axis
  scaling BENCH_JSON     CPU single-thread, CPU multi-thread and CUDA time per iteration
  wall RUN_DIR           wall shear and heat-flux distributions from wall.csv
  sweep SWEEP_JSON       2D thrust against the quasi-1D isentropic estimate
Every figure states grid, backend/precision and termination where applicable.
"""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import xml.etree.ElementTree as ET
os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[1] / ".cache/matplotlib"))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# Validated reference palette (dataviz skill): categorical slots and the one-hue blue ramp.
BLUE, ORANGE, AQUA, YELLOW = '#2a78d6', '#eb6834', '#1baf7a', '#eda100'
INK, INK2, MUTED, GRID, AXIS, SURFACE = '#0b0b0b', '#52514e', '#898781', '#e1e0d9', '#c3c2b7', '#fcfcfb'
RAMP = ['#cde2fb', '#b7d3f6', '#9ec5f4', '#86b6ef', '#6da7ec', '#5598e7', '#3987e5', '#2a78d6',
        '#256abf', '#1c5cab', '#184f95', '#104281', '#0d366b']
SEQUENTIAL = LinearSegmentedColormap.from_list('astraflow_blue', RAMP)


def style():
    plt.rcParams.update({
        'font.family': 'sans-serif', 'font.size': 10, 'text.color': INK, 'axes.labelcolor': INK2,
        'axes.edgecolor': AXIS, 'axes.linewidth': 1, 'axes.spines.top': False, 'axes.spines.right': False,
        'xtick.color': MUTED, 'ytick.color': MUTED, 'grid.color': GRID, 'grid.linewidth': 1,
        'axes.grid': True, 'axes.axisbelow': True, 'lines.linewidth': 2, 'lines.solid_joinstyle': 'round',
        'lines.solid_capstyle': 'round', 'legend.frameon': False, 'figure.facecolor': SURFACE,
        'axes.facecolor': SURFACE, 'savefig.dpi': 180, 'savefig.facecolor': SURFACE,
        'figure.constrained_layout.use': True})


def read_vts(path):
    root = ET.parse(path).getroot()
    piece = root.find('./StructuredGrid/Piece')
    nx, nr = (int(v) for v in piece.attrib['Extent'].split()[1:4:2])
    pts = [float(v) for v in piece.find('./Points/DataArray').text.split()]
    x = [pts[3 * i] for i in range((nx + 1) * (nr + 1))]
    r = [pts[3 * i + 1] for i in range((nx + 1) * (nr + 1))]
    arrays = {a.attrib['Name']: [float(v) for v in a.text.split()] for a in piece.findall('./CellData/DataArray')}
    return nx, nr, x, r, arrays


def exact_sod(xs, t, gamma=1.4, x0=0.5):
    """Independent exact Riemann solution (Toro): star pressure by Newton on the pressure function."""
    rl, ul, pl, rr, ur, pr = 1.0, 0.0, 1.0, 0.125, 0.0, 0.1
    al, ar = math.sqrt(gamma * pl / rl), math.sqrt(gamma * pr / rr)
    g1, g2, g3, g4 = (gamma - 1) / (2 * gamma), (gamma + 1) / (2 * gamma), 2 / (gamma - 1), 2 / (gamma + 1)

    def f(p, rk, pk, ak):
        if p > pk:
            A, B = g4 / rk, (gamma - 1) / (gamma + 1) * pk
            return (p - pk) * math.sqrt(A / (p + B)), math.sqrt(A / (B + p)) * (1 - (p - pk) / (2 * (B + p)))
        return g3 * ak * ((p / pk) ** g1 - 1), (1 / (rk * ak)) * (p / pk) ** (-g2)
    p = 0.5 * (pl + pr)
    for _ in range(100):
        fl, dl = f(p, rl, pl, al)
        fr, dr = f(p, rr, pr, ar)
        step = (fl + fr + ur - ul) / (dl + dr)
        p = max(p - step, 1e-12)
        if abs(step) < 1e-14:
            break
    ps = p
    us = 0.5 * (ul + ur) + 0.5 * (f(ps, rr, pr, ar)[0] - f(ps, rl, pl, al)[0])
    out = []
    for x in xs:
        s = (x - x0) / t
        if s <= us:  # left of contact
            if ps > pl:  # left shock
                sl = ul - al * math.sqrt(g2 * ps / pl + g1)
                if s < sl:
                    out.append((rl, ul, pl))
                else:
                    out.append((rl * (ps / pl + (gamma - 1) / (gamma + 1)) / ((gamma - 1) / (gamma + 1) * ps / pl + 1), us, ps))
            else:  # left rarefaction
                shl, asl = ul - al, al * (ps / pl) ** g1
                stl = us - asl
                if s < shl:
                    out.append((rl, ul, pl))
                elif s > stl:
                    out.append((rl * (ps / pl) ** (1 / gamma), us, ps))
                else:
                    a = g4 * (al + (gamma - 1) / 2 * (ul - s))
                    out.append((rl * (a / al) ** g3, g4 * (al + (gamma - 1) / 2 * ul + s), pl * (a / al) ** (1 / g1)))
        else:
            if ps > pr:  # right shock
                sr = ur + ar * math.sqrt(g2 * ps / pr + g1)
                if s > sr:
                    out.append((rr, ur, pr))
                else:
                    out.append((rr * (ps / pr + (gamma - 1) / (gamma + 1)) / ((gamma - 1) / (gamma + 1) * ps / pr + 1), us, ps))
            else:
                shr, asr = ur + ar, ar * (ps / pr) ** g1
                if s > shr:
                    out.append((rr, ur, pr))
                elif s < us + asr:
                    out.append((rr * (ps / pr) ** (1 / gamma), us, ps))
                else:
                    a = g4 * (ar - (gamma - 1) / 2 * (ur - s))
                    out.append((rr * (a / ar) ** g3, g4 * (-ar + (gamma - 1) / 2 * ur + s), pr * (a / ar) ** (1 / g1)))
    return out, ps, us


def sod(run, out):
    nx, nr, x, r, arrays = read_vts(run / 'final_state.vts')
    summary = json.loads((run / 'summary.json').read_text())
    t = summary['simulated_time']
    xc = [0.5 * (x[i] + x[i + 1]) for i in range(nx)]
    rho = arrays['density'][:nx]
    u = [arrays['velocity'][3 * i] for i in range(nx)]
    p = arrays['pressure'][:nx]
    fine = [i / 2000 for i in range(2001)]
    exact, ps, us = exact_sod(fine, t)
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.4))
    errors = []
    for ax, k, name, num in zip(axes, range(3), ('Density', 'Velocity', 'Pressure'), (rho, u, p)):
        ax.plot(fine, [e[k] for e in exact], color=INK2, lw=1.2, label='Exact Riemann solution')
        ax.plot(xc, num, color=BLUE, lw=0, marker='o', ms=2.6, label=f'AstraFlow, {nx} cells')
        ax.set_xlabel('x')
        ax.set_title(name, loc='left', color=INK)
        ax.set_xlim(0, 1)
        ex, _, _ = exact_sod(xc, t)
        errors.append(sum(abs(a - e[k]) for a, e in zip(num, ex)) / nx)
    axes[0].legend(loc='upper right', fontsize=8)
    fig.suptitle(f'Sod shock tube at t = {t:g}, {summary["backend"]} {summary["precision"]}, MUSCL-MC + HLLC + SSP-RK2'
                 f'   (L1 errors: rho {errors[0]:.2e}, u {errors[1]:.2e}, p {errors[2]:.2e}; p* = {ps:.5f})',
                 fontsize=9, color=INK2, x=0.01, ha='left')
    fig.savefig(out / 'sod-verification.png')
    plt.close(fig)
    return {'L1': dict(zip(('rho', 'u', 'p'), errors)), 'p_star': ps, 'u_star': us}


def field(run, out):
    nx, nr, x, r, arrays = read_vts(run / 'final_state.vts')
    summary = json.loads((run / 'summary.json').read_text())
    import numpy as np
    X = np.array(x).reshape(nr + 1, nx + 1)
    R = np.array(r).reshape(nr + 1, nx + 1)
    label = f'{nx} × {nr}, {summary["backend"]} {summary["precision"]}, {summary["termination_reason"]}'
    for name, key, unit, scale in (('nozzle-mach', 'Mach', 'Mach number', 1), ('nozzle-pressure', 'pressure', 'Pressure [MPa]', 1e-6),
                                   ('nozzle-temperature', 'temperature', 'Temperature [K]', 1)):
        C = np.array(arrays[key]).reshape(nr, nx) * scale
        fig, ax = plt.subplots(figsize=(10, 3.4))
        magnification = 4
        for sign in (1, -1):
            m = ax.pcolormesh(X * 1e3, sign * R * 1e3 * magnification, C, cmap=SEQUENTIAL, shading='flat', rasterized=True)
        ax.plot(X[-1] * 1e3, R[-1] * 1e3 * magnification, color=INK, lw=1)
        ax.plot(X[-1] * 1e3, -R[-1] * 1e3 * magnification, color=INK, lw=1)
        ax.set_aspect('equal')
        ax.grid(False)
        ax.set_xlabel('x [mm]')
        ax.set_ylabel(f'r [mm], radial axis magnified {magnification}x')
        ticks = [-6, -3, 0, 3, 6]
        ax.set_yticks([t * magnification for t in ticks])
        ax.set_yticklabels([str(t) for t in ticks])
        ax.set_title(f'{unit} — {label}', loc='left', color=INK)
        cb = fig.colorbar(m, ax=ax, orientation='horizontal', shrink=0.45, pad=0.22, aspect=40)
        cb.outline.set_visible(False)
        fig.savefig(out / f'{name}.png')
        plt.close(fig)


def scaling(bench, out):
    data = json.loads(Path(bench).read_text())
    rows = data['results']
    cells = [r['cells'] for r in rows]
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    series = [('CPU, 1 thread', [r['cpu_iteration_ms'] for r in rows], BLUE),
              (f'CPU, {rows[0].get("cpu_threads", "n")} threads', [r.get('cpu_parallel_iteration_ms') for r in rows], ORANGE),
              ('CUDA FP32, RTX 5070 Laptop', [r.get('gpu_iteration_ms') for r in rows], AQUA)]
    for name, values, color in series:
        if any(v is None for v in values):
            continue
        ax.plot(cells, values, color=color, marker='o', ms=5, label=name,
                markeredgecolor=SURFACE, markeredgewidth=1.5)
        ax.annotate(f'{values[-1]:.2f} ms', (cells[-1], values[-1]), xytext=(6, 0), textcoords='offset points',
                    va='center', fontsize=8, color=INK2)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Cells')
    ax.set_ylabel('Time per iteration [ms]')
    ax.set_xticks(cells)
    ax.set_xticklabels([f'{r["nx"]}×{r["nr"]}' for r in rows])
    ax.set_title('Viscous nozzle, equal FP32, 20 warm-up + measured iterations; sample means, not guaranteed speedups',
                 loc='left', fontsize=9, color=INK2)
    ax.legend(loc='upper left')
    ax.set_xlim(cells[0] * 0.8, cells[-1] * 1.9)
    fig.savefig(out / 'cpu-gpu-scaling.png')
    plt.close(fig)


def wall(run, out):
    rows = list(csv.DictReader((run / 'wall.csv').open()))
    summary = json.loads((run / 'summary.json').read_text())
    config = json.loads((run / 'config.json').read_text())
    L = sum(config['geometry'][k] for k in ('chamber_length', 'contraction_length', 'expansion_length'))
    x = [float(r['x']) / L for r in rows]
    fig, axes = plt.subplots(2, 1, figsize=(8, 5.2), sharex=True)
    axes[0].plot(x, [-float(r['traction_x']) for r in rows], color=BLUE)
    axes[0].set_ylabel('Wall shear stress [Pa]')
    axes[1].plot(x, [float(r['heat_flux_into_wall']) / 1e6 for r in rows], color=ORANGE)
    axes[1].set_ylabel('Heat flux into wall [MW/m²]')
    axes[1].set_xlabel('Axial position, x/L')
    throat = (config['geometry']['chamber_length'] + config['geometry']['contraction_length']) / L
    for ax in axes:
        ax.axvline(throat, color=AXIS, lw=1)
    axes[0].annotate('throat', (throat, axes[0].get_ylim()[1]), xytext=(4, -10), textcoords='offset points', fontsize=8, color=MUTED)
    w = summary.get('wall', {})
    axes[0].set_title(f'Cooled nozzle wall, T_wall = {config["boundary_conditions"]["wall_temperature"]:g} K, Sutherland viscosity — '
                      f'{summary["grid"][0]} × {summary["grid"][1]}, {summary["backend"]} {summary["precision"]}, {summary["termination_reason"]}\n'
                      f'Total wall heat transfer {w.get("wall_heat_transfer_rate", float("nan")) / 1e3:.2f} kW, '
                      f'viscous drag {w.get("wall_axial_viscous_force", float("nan")):.3f} N', loc='left', fontsize=9, color=INK2)
    axes[0].set_xlim(0, 1)
    fig.savefig(out / 'wall-heat-flux.png')
    plt.close(fig)


def sweep(path, out):
    data = json.loads(Path(path).read_text())
    cases = [c for c in data['cases'] if c.get('converged')]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    ar = [c['area_ratio'] for c in cases]
    axes[0].plot(ar, [c['thrust_1d'] for c in cases], color=INK2, lw=1.2, label='Quasi-1D isentropic')
    axes[0].plot(ar, [c['thrust_2d'] for c in cases], color=BLUE, lw=0, marker='o', ms=6, label='AstraFlow 2D axisymmetric',
                 markeredgecolor=SURFACE, markeredgewidth=1.5)
    axes[0].set_xlabel('Exit-to-throat area ratio')
    axes[0].set_ylabel('Thrust [N]')
    axes[0].legend()
    axes[1].plot(ar, [100 * c['thrust_relative_difference'] for c in cases], color=BLUE, marker='o', ms=5, markeredgecolor=SURFACE, markeredgewidth=1.5)
    axes[1].axhline(0, color=AXIS, lw=1)
    axes[1].set_xlabel('Exit-to-throat area ratio')
    axes[1].set_ylabel('2D minus quasi-1D thrust [%]')
    fig.suptitle(f'Viscous no-slip nozzle expansion sweep against inviscid quasi-1D theory, {len(cases)} converged steady cases (local pseudo-time, van Albada, CUDA FP64)',
                 fontsize=9, color=INK2, x=0.01, ha='left')
    fig.savefig(out / 'nozzle-sweep.png')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('kind', choices=('sod', 'field', 'scaling', 'wall', 'sweep'))
    parser.add_argument('source', type=Path)
    parser.add_argument('--output', type=Path, default=Path(__file__).resolve().parents[1] / 'docs/figures')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    style()
    result = {'sod': sod, 'field': field, 'scaling': scaling, 'wall': wall, 'sweep': sweep}[args.kind](args.source, args.output)
    if result:
        print(json.dumps(result, indent=2))
    print(f'Wrote {args.kind} figure(s) to {args.output}')


if __name__ == '__main__':
    main()

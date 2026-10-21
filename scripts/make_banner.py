#!/usr/bin/env python3
"""Render the README banner as a standalone SVG from the actual example nozzle contour."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def radius(g, x):
    """Same C1 smoothstep contour as geometry/nozzle.cpp."""
    if x <= g['chamber_length']:
        return g['chamber_radius']
    contraction = x < g['chamber_length'] + g['contraction_length']
    if contraction:
        t = (x - g['chamber_length']) / g['contraction_length']
        a, b = g['chamber_radius'], g['throat_radius']
    else:
        t = (x - g['chamber_length'] - g['contraction_length']) / g['expansion_length']
        a, b = g['throat_radius'], g['exit_radius']
    t = min(max(t, 0.0), 1.0)
    return a + (b - a) * t * t * (3 - 2 * t)


def main():
    g = json.loads((ROOT / 'examples/rocket_nozzle/config.json').read_text())['geometry']
    L = g['chamber_length'] + g['contraction_length'] + g['expansion_length']
    W, H = 1200, 320
    x0, x1, cy = 640, 1150, H / 2
    scale_x = (x1 - x0) / L
    scale_r = 150 / g['exit_radius']
    n = 160
    upper = [(x0 + scale_x * (L * i / n), cy - scale_r * radius(g, L * i / n)) for i in range(n + 1)]
    lower = [(x, 2 * cy - y) for x, y in reversed(upper)]
    path = 'M ' + ' L '.join(f'{x:.1f},{y:.1f}' for x, y in upper + lower) + ' Z'
    # Mach-like shading bands: stations coloured by the quasi-1D area ratio, light to dark.
    bands = []
    for i in range(n):
        xa, xb = L * i / n, L * (i + 1) / n
        ra, rb = radius(g, xa), radius(g, xb)
        sup = xb > g['chamber_length'] + g['contraction_length']
        ratio = (0.5 * (ra + rb) / g['throat_radius']) ** 2
        level = min(1.0, 0.15 + 0.85 * (0.35 * (1 - 1 / ratio) if not sup else 0.45 + 0.55 * (ratio - 1) / ((g['exit_radius'] / g['throat_radius']) ** 2 - 1)))
        bands.append((x0 + scale_x * xa, x0 + scale_x * xb, cy - scale_r * ra, cy - scale_r * rb, level))
    rects = []
    for xa, xb, ya, yb, level in bands:
        top = min(ya, yb)
        rects.append(f'<rect x="{xa:.1f}" y="{top:.1f}" width="{xb - xa + 0.6:.1f}" height="{2 * (cy - top):.1f}" fill="#9ec5f4" fill-opacity="{0.25 + 0.7 * level:.2f}"/>')
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" role="img" aria-label="AstraFlow banner">
  <defs>
    <linearGradient id="bg" x1="0" y1="0" x2="1" y2="1">
      <stop offset="0" stop-color="#0b1220"/>
      <stop offset="1" stop-color="#10263f"/>
    </linearGradient>
    <clipPath id="nozzle"><path d="{path}"/></clipPath>
    <linearGradient id="glow" x1="0" y1="0" x2="1" y2="0">
      <stop offset="0" stop-color="#2a78d6" stop-opacity="0.0"/>
      <stop offset="1" stop-color="#2a78d6" stop-opacity="0.35"/>
    </linearGradient>
  </defs>
  <rect width="{W}" height="{H}" fill="url(#bg)"/>
  <g clip-path="url(#nozzle)">
    <rect x="{x0}" y="0" width="{x1 - x0}" height="{H}" fill="#0d366b"/>
    {''.join(rects)}
    <rect x="{x0}" y="0" width="{x1 - x0}" height="{H}" fill="url(#glow)"/>
  </g>
  <path d="{path}" fill="none" stroke="#cde2fb" stroke-width="2"/>
  <line x1="{x0}" y1="{cy}" x2="{x1}" y2="{cy}" stroke="#86b6ef" stroke-width="1" stroke-dasharray="6 6" opacity="0.6"/>
  <text x="60" y="128" font-family="system-ui, -apple-system, 'Segoe UI', sans-serif" font-size="64" font-weight="700" fill="#f0efec">AstraFlow</text>
  <text x="62" y="170" font-family="system-ui, -apple-system, 'Segoe UI', sans-serif" font-size="21" fill="#c3c2b7">GPU-accelerated axisymmetric compressible CFD</text>
  <text x="62" y="200" font-family="system-ui, -apple-system, 'Segoe UI', sans-serif" font-size="21" fill="#c3c2b7">for rocket nozzle flows</text>
  <g font-family="system-ui, -apple-system, 'Segoe UI', sans-serif" font-size="14" fill="#cde2fb">
    <rect x="62" y="232" width="92" height="30" rx="15" fill="#1c5cab" fill-opacity="0.55"/><text x="108" y="252" text-anchor="middle">C++20</text>
    <rect x="164" y="232" width="128" height="30" rx="15" fill="#1c5cab" fill-opacity="0.55"/><text x="228" y="252" text-anchor="middle">CUDA sm_120</text>
    <rect x="302" y="232" width="150" height="30" rx="15" fill="#1c5cab" fill-opacity="0.55"/><text x="377" y="252" text-anchor="middle">Verified numerics</text>
    <rect x="462" y="232" width="118" height="30" rx="15" fill="#1c5cab" fill-opacity="0.55"/><text x="521" y="252" text-anchor="middle">MIT licence</text>
  </g>
  <text x="{x1}" y="{H - 14}" text-anchor="end" font-family="system-ui, -apple-system, 'Segoe UI', sans-serif" font-size="12" fill="#898781">Example nozzle contour: chamber 6 mm, throat 4.5 mm, exit 6.75 mm, length 210 mm</text>
</svg>
'''
    out = ROOT / 'docs/assets/banner.svg'
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(svg)
    print(f'Wrote {out}')
    pipeline()


def pipeline():
    """Solver pipeline diagram: one timestep of the shared CPU/CUDA finite-volume operator."""
    stages = [('Conservative state', 'rho, rho u, rho v, rho E (SoA)'),
              ('Primitive conversion', 'positivity check, CFL bound per cell'),
              ('Reconstruction', 'MUSCL slopes + transport gradients'),
              ('Face fluxes', 'rotated HLLC / HLLE + corrected viscous flux'),
              ('Residual', 'face sums + axisymmetric source'),
              ('SSP-RK2 update', 'two stages, or local pseudo-time'),
              ('Diagnostics', 'residual RMS, boundary mass flux, one transfer')]
    W, H, bw, bh, gap = 1200, 190, 150, 92, 22
    x0 = (W - (len(stages) * bw + (len(stages) - 1) * gap)) / 2
    parts = ['<defs><linearGradient id="g" x1="0" y1="0" x2="1" y2="1"><stop offset="0" stop-color="#0b1220"/><stop offset="1" stop-color="#10263f"/></linearGradient>'
             '<marker id="a" markerWidth="8" markerHeight="8" refX="6" refY="4" orient="auto"><path d="M0,0 L8,4 L0,8 z" fill="#86b6ef"/></marker></defs>',
             f'<rect width="{W}" height="{H}" fill="url(#g)"/>']
    for i, (title, sub) in enumerate(stages):
        x = x0 + i * (bw + gap)
        y = (H - bh) / 2
        parts.append(f'<rect x="{x:.0f}" y="{y:.0f}" width="{bw}" height="{bh}" rx="12" fill="#1c5cab" fill-opacity="0.45" stroke="#86b6ef" stroke-width="1.2"/>')
        parts.append(f'<text x="{x + bw / 2:.0f}" y="{y + 34:.0f}" text-anchor="middle" font-family="system-ui, -apple-system, \'Segoe UI\', sans-serif" font-size="14" font-weight="600" fill="#f0efec">{title}</text>')
        words, lines, cur = sub.split(), [], ''
        for w in words:
            if len(cur) + len(w) + 1 > 22:
                lines.append(cur)
                cur = w
            else:
                cur = (cur + ' ' + w).strip()
        lines.append(cur)
        for j, line in enumerate(lines[:3]):
            parts.append(f'<text x="{x + bw / 2:.0f}" y="{y + 54 + 15 * j:.0f}" text-anchor="middle" font-family="system-ui, -apple-system, \'Segoe UI\', sans-serif" font-size="11" fill="#c3c2b7">{line}</text>')
        if i < len(stages) - 1:
            parts.append(f'<line x1="{x + bw + 3:.0f}" y1="{H / 2:.0f}" x2="{x + bw + gap - 3:.0f}" y2="{H / 2:.0f}" stroke="#86b6ef" stroke-width="1.5" marker-end="url(#a)"/>')
    parts.append(f'<text x="{W / 2:.0f}" y="{H - 12}" text-anchor="middle" font-family="system-ui, -apple-system, \'Segoe UI\', sans-serif" font-size="12" fill="#898781">One timestep of the shared header-only operator; CPU and CUDA call the same functions and differ only in scheduling and memory</text>')
    svg = f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" role="img" aria-label="AstraFlow solver pipeline">' + ''.join(parts) + '</svg>\n'
    out = ROOT / 'docs/assets/pipeline.svg'
    out.write_text(svg)
    print(f'Wrote {out}')


if __name__ == '__main__':
    main()

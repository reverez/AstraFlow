#!/usr/bin/env python3
"""Validate an exported AstraFlow run using only the Python standard library."""
import csv
import json
import math
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

def validate(directory):
    root = Path(directory)
    config = json.loads((root / 'config.json').read_text())
    summary = json.loads((root / 'summary.json').read_text())
    performance = json.loads((root / 'performance.json').read_text())
    nx, nr = summary['grid']
    assert [config['mesh']['nx'], config['mesh']['nr']] == [nx, nr]
    assert summary['iterations'] > 0 and performance['wall_ms'] >= 0
    for filename in ('mesh.vts', 'final_state.vts'):
        tree = ET.parse(root / filename)
        assert tree.getroot().attrib['type'] == 'StructuredGrid'
        piece = tree.find('./StructuredGrid/Piece')
        assert piece.attrib['Extent'] == f'0 {nx} 0 {nr} 0 0'
        values = piece.find('./Points/DataArray').text.split()
        assert len(values) == 3 * (nx + 1) * (nr + 1)
        arrays = piece.findall('./CellData/DataArray')
        if filename == 'final_state.vts':
            assert {a.attrib['Name'] for a in arrays} == {'density','pressure','temperature','velocity','Mach','total_energy'}
        for a in arrays:
            values = [float(x) for x in a.text.split()]
            assert len(values) == nx * nr * int(a.attrib['NumberOfComponents'])
            assert all(math.isfinite(x) for x in values)
            if a.attrib['Name'] in ('density','pressure','temperature'):
                assert min(values) > 0
    rows = list(csv.DictReader((root / 'convergence.csv').open()))
    assert int(rows[-1]['iteration']) == summary['iterations']
    if summary['specific_impulse'] is not None:
        assert math.isclose(summary['specific_impulse'] * summary['mass_flow'] * 9.80665, summary['estimated_thrust'], rel_tol=1e-12)
    print(f'{root}: JSON/CSV/VTK structure, finite positive fields and engineering consistency PASS')

if __name__ == '__main__':
    for directory in sys.argv[1:]:
        validate(directory)

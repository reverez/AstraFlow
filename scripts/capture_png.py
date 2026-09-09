#!/usr/bin/env python3
"""Losslessly encode AstraFlow's actual PPM framebuffer capture as PNG; no dependencies."""
from pathlib import Path
import struct
import sys
import zlib

def convert(source, destination):
    with Path(source).open('rb') as f:
        assert f.readline().strip() == b'P6'
        width, height = map(int, f.readline().split())
        assert f.readline().strip() == b'255'
        pixels = f.read()
    assert len(pixels) == width * height * 3
    assert len(set(pixels[i:i+3] for i in range(0, len(pixels), 3))) > 100, 'Capture appears blank'
    def chunk(kind, data):
        return struct.pack('>I', len(data)) + kind + data + struct.pack('>I', zlib.crc32(kind + data))
    rows = b''.join(b'\x00' + pixels[y*width*3:(y+1)*width*3] for y in range(height))
    png = b'\x89PNG\r\n\x1a\n' + chunk(b'IHDR', struct.pack('>IIBBBBB', width, height, 8, 2, 0, 0, 0))
    png += chunk(b'IDAT', zlib.compress(rows, 9)) + chunk(b'IEND', b'')
    target = Path(destination)
    if target.exists():
        raise FileExistsError(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(png)
    print(f'Encoded actual {width}x{height} framebuffer: {target}')

if __name__ == '__main__':
    convert(sys.argv[1], sys.argv[2])

#!/usr/bin/env python3
"""Install only CUDA 12.8.1 compiler/runtime/CCCL into this project's .toolchains.

Uses NVIDIA's pinned SHA256 hashes; needs network access, never sudo.
CUDA is distributed under NVIDIA's license, not AstraFlow's MIT license.
"""
import hashlib
from pathlib import Path
import tarfile
import urllib.request

ROOT = Path(__file__).resolve().parents[1] / '.toolchains'
PACKAGES = {
    'cuda_nvcc/linux-x86_64/cuda_nvcc-linux-x86_64-12.8.93-archive.tar.xz':
        '9961b3484b6b71314063709a4f9529654f96782ad39e72bf1e00f070db8210d3',
    'cuda_cudart/linux-x86_64/cuda_cudart-linux-x86_64-12.8.90-archive.tar.xz':
        '8d566b5fe745c46842dc16945cf36686227536decd2302c372be86da37faca68',
    'cuda_cccl/linux-x86_64/cuda_cccl-linux-x86_64-12.8.90-archive.tar.xz':
        '0740e9e01e4f15e17c5ab8d68bba4f8ec0eb6b84edccba4ac45112d2d2174e4b',
}

def main():
    downloads = ROOT / 'downloads'
    downloads.mkdir(parents=True, exist_ok=True)
    destination = ROOT / 'cuda-12.8.1'
    destination.mkdir(exist_ok=True)
    for name, digest in PACKAGES.items():
        archive = downloads / Path(name).name
        if not archive.exists():
            urllib.request.urlretrieve('https://developer.download.nvidia.com/compute/cuda/redist/' + name, archive)
        if hashlib.sha256(archive.read_bytes()).hexdigest() != digest:
            raise RuntimeError(f'CUDA archive checksum mismatch: {archive.name}')
        with tarfile.open(archive) as source:
            for member in source.getmembers():
                parts = Path(member.name).parts
                if len(parts) > 1:
                    member.name = str(Path(*parts[1:]))
                    source.extract(member, destination, filter='data')
        print(f'Installed and verified {archive.name}')
    if not (destination / 'lib64').exists():
        (destination / 'lib64').symlink_to('lib', target_is_directory=True)
    print(f'CMake CUDA compiler: {destination / "bin/nvcc"}')

if __name__ == '__main__':
    main()

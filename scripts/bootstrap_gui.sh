#!/usr/bin/env bash
# Ubuntu 24.04 project-local GUI headers. No sudo or system package installation.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p .toolchains/gui-packages .toolchains/gui-sysroot
cd .toolchains/gui-packages
apt-get download \
  libgl-dev=1.7.0-1build1 libglx-dev=1.7.0-1build1 libglvnd-core-dev=1.7.0-1build1 \
  libx11-dev=2:1.8.7-1build1 libxrandr-dev=2:1.5.2-2build1 \
  libxinerama-dev=2:1.1.4-3build1 libxcursor-dev=1:1.2.1-1build1 \
  libxi-dev=2:1.8.1-1build1 libxext-dev=2:1.3.4-1build2 \
  libxfixes-dev=1:6.0.0-2build1 libxrender-dev=1:0.9.10-1.1build1 \
  x11proto-dev=2023.2-1 libxau-dev=1:1.0.9-1build6 \
  libxdmcp-dev=1:1.1.3-0ubuntu6 libxcb1-dev=1.15-1ubuntu2 xtrans-dev=1.4.0-1
for package in *.deb; do dpkg-deb -x "$package" ../gui-sysroot; done
python3 - <<'PY'
from pathlib import Path
root = Path('../gui-sysroot/usr/lib/x86_64-linux-gnu')
for library in root.glob('*.so'):
    if library.is_symlink() and not library.exists():
        target = library.readlink()
        system = Path('/usr/lib/x86_64-linux-gnu') / target
        if len(target.parts) == 1 and system.exists():
            (root / target).symlink_to(system)
print('GUI development files installed locally; existing system runtime libraries are reused.')
PY

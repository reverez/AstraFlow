# Development environment

Preflight performed once on 2026-09-09 at `/home/arnav/dev/projects/portfolio`; canonical project created at `/home/arnav/dev/projects/AstraFlow`. No existing project data was present.

- WSL2 kernel: 6.18.33.2-microsoft-standard-WSL2, x86_64.
- Ubuntu 24.04.4 LTS.
- GCC 13.3.0, CMake 3.28.3, Ninja 1.11.1, Git 2.43.0, Python 3.12.3.
- Initial sandbox NVML call: GPU access blocked by the operating system.
- One outside-sandbox diagnostic succeeded: NVIDIA GeForce RTX 5070 family, 8151 MiB, driver 615.65.06, host KMD 616.56, driver-reported CUDA UMD 13.4. Driver capability is not an installed toolkit version.
- No installed nvcc at preflight. Ubuntu package candidate CUDA 12.0 does not meet SM120 requirements. sudo needs a password.
- User approved project-local dependency/toolkit downloads. CUDA 12.8.1 compiler/runtime/CCCL redistribution packages are SHA256-checked against NVIDIA's manifest and installed under ignored `.toolchains/cuda-12.8.1`; no system driver changes.
- GitHub CLI initially absent; gh 2.65.0 installed project-locally, but no host is authenticated.
- DISPLAY=:0 and WAYLAND_DISPLAY=wayland-0; actual WSLg GUI startup, rendering and scripted interaction subsequently passed.
- clang-format is available; nsys, ncu and compute-sanitizer were not found on PATH.

Sources: [NVIDIA CUDA 12.8.1 manifest](https://developer.download.nvidia.com/compute/cuda/redist/redistrib_12.8.1.json), [NVIDIA Linux installation guide](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html).

Milestone 0 runtime probe: NVIDIA GeForce RTX 5070 Laptop GPU; compute capability 12.0; global memory 8,518,041,600 bytes; CUDA runtime 12080; executed kernel architecture sm_120. NVCC 12.8.93. CMake 3.28.3 accepts architecture 120 explicitly on the target. Local gh 2.65.0 reports no authenticated hosts.

Relocation update: the user moved the repository to `/home/arnav/dev/projects/portfolio/AstraFlow`. Continued development there on 2026-09-09. The old absolute-path CMake trees are preserved as `build-before-move` and `build-cpu-before-move`; fresh build trees reuse their downloaded dependency sources. No project data or user commits were removed.

Benchmark CPU: Intel Core Ultra 9 285H, 16 logical CPUs visible to WSL; CPU reference solver uses one thread. A project `.venv` uses WSL Python 3.12.3; auxiliary scripts need no installed packages.

Final remote update: origin is `git@github.com:reverez/AstraFlow.git`; SSH fetch succeeds. Four user README commits were discovered and merged normally before final publication. The earlier unauthenticated gh status did not prevent SSH Git access.

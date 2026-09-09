# AstraFlow

AstraFlow is being implemented as a GPU-accelerated axisymmetric compressible CFD simulator for rocket nozzle flows, using C++20 and CUDA, with a CPU reference backend.

Development Milestones are underway

**Development status:** Milestones 0–6 complete: verified CPU/CUDA Euler and axisymmetric viscous finite volumes, CLI, engineering analysis and reproducible output. GUI/performance/final review remain in progress. See [Phase 1 status](docs/PHASE1_STATUS.md).

## Build

Dependencies are pinned: nlohmann/json 3.11.3 and Catch2 3.7.1 (MIT and BSL-1.0 respectively; upstream licenses remain in fetched sources). CMake downloads them at first configure. CPU builds need CMake 3.24+, Ninja and a C++20 compiler.


The initial device check executes a CUDA kernel and reports its compiled architecture. WSL GPU access may require execution outside an agent sandbox.

## Scientific scope

The planned V1 model is a single-species, calorically perfect ideal gas with conservative finite volumes, MUSCL, HLLC/HLLE, SSP-RK2, axisymmetric geometry and Newtonian viscosity/Fourier heat conduction. This is a numerical discretisation of the Euler/Navier–Stokes equations; it is not a solution of the mathematical existence/smoothness problem.

Combustion, turbulence, external plumes and advanced thermodynamics are outside V1.

MIT license. CUDA remains subject to NVIDIA's license.

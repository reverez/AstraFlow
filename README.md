# AstraFlow

> **GPU-accelerated compressible CFD for axisymmetric rocket nozzle flows.**

AstraFlow is a scientific-computing project implementing a **2D axisymmetric compressible Euler/Navier–Stokes solver** in **C++20 and NVIDIA CUDA**, with a CPU reference backend, scientific verification suite, engineering analysis pipeline, reproducible simulation output, and interactive visualisation.

The project targets an **NVIDIA GeForce RTX 5070 Laptop GPU** using NVIDIA Blackwell **compute capability 12.0 (`sm_120`)**.

I architected the mathematical formulation, numerical method, physical assumptions, verification strategy, CPU/GPU architecture, project structure, and milestone acceptance criteria. AI coding agents are used to accelerate implementation, testing, debugging, profiling, and documentation under those specifications.

Numerical results are accepted only after comparison against analytical solutions, conservation properties, convergence behaviour, or CPU/GPU parity tests.

---

## Status

**Phase 1 / V1 implementation and local verification are complete.**

Milestones **0–9 are complete**, including verified CPU/CUDA Euler and axisymmetric Navier–Stokes solvers, engineering analysis, reproducible CLI output, the interactive GUI, CUDA profiling/optimisation, grid refinement and final acceptance tests.

Final checks passed: Release CPU **6/6**, Release CUDA **8/8**, CPU AddressSanitizer/UndefinedBehaviorSanitizer with leak detection **6/6**, both 3,000-step rocket workflows, independent output validation and actual WSLg GUI smoke/capture. The project is ready for architectural review.

See [`docs/PHASE1_STATUS.md`](docs/PHASE1_STATUS.md).

---

## Preview

<!-- Further standalone field and scaling figures may be added later. -->

### Interactive Simulation

![Actual AstraFlow interface on WSLg](docs/screenshots/astraflow.png)

This is an actual framebuffer capture with nozzle flow, live residual history, engineering results and CUDA telemetry. Radial display magnification is labelled and does not modify physical geometry.

### Mach Field

<!--
<p align="center">
  <img src="docs/assets/nozzle-mach.png" width="95%" alt="AstraFlow Mach-number field">
</p>
-->

### Pressure Field

<!--
<p align="center">
  <img src="docs/assets/nozzle-pressure.png" width="95%" alt="AstraFlow pressure field">
</p>
-->

### CPU vs CUDA Scaling

<!--
<p align="center">
  <img src="docs/assets/cpu-gpu-scaling.png" width="90%" alt="CPU versus CUDA scaling">
</p>
-->

---

# Objectives

AstraFlow investigates the complete path from a continuous physical model to verified high-performance numerical software.

The principal objectives are to:

1. implement the governing equations rather than wrap an existing CFD package;
2. construct a conservative numerical method for compressible and transonic flow;
3. reproduce characteristic internal rocket-nozzle behaviour;
4. verify the implementation quantitatively;
5. maintain CPU and CUDA backends for cross-validation;
6. exploit modern GPU hardware for parallel finite-volume computation;
7. measure numerical error and hardware scaling;
8. expose the solver through reproducible CLI and interactive interfaces.

---

# Authorship and AI-Assisted Implementation

AstraFlow follows an **architect-directed, AI-assisted engineering workflow**.

I am responsible for the design and technical review of:

- physical scope;
- governing equations;
- conservative state representation;
- finite-volume discretisation;
- reconstruction strategy;
- Riemann solver;
- temporal integration;
- axisymmetric formulation;
- viscous and thermal transport;
- boundary conditions;
- numerical stability;
- verification cases;
- error tolerances;
- CPU/CUDA parity strategy;
- CUDA architecture;
- software architecture;
- milestone gates;
- interpretation of numerical results.

AI coding agents, principally Codex, are used to accelerate:

- C++ implementation;
- CUDA implementation;
- unit and verification tests;
- build-system development;
- debugging;
- targeted numerical experiments;
- profiling;
- repetitive infrastructure;
- technical documentation.

The workflow is therefore:

```text
Mathematical architecture
        ↓
Implementation specification
        ↓
AI-assisted implementation
        ↓
Numerical verification
        ↓
Architectural review
        ↓
Milestone acceptance
```

Generated code is not accepted solely because it compiles or produces plausible flow fields.

---

# Physical Model

V1 models a **single-species, calorically perfect ideal gas** flowing through an axisymmetric converging-diverging nozzle.

The conservative state is

```math
\mathbf{U}
=
\begin{bmatrix}
\rho \\
\rho u \\
\rho v \\
\rho E
\end{bmatrix}
```

where:

- $\rho$ — density;
- $u$ — axial velocity;
- $v$ — radial velocity;
- $E$ — total specific energy.

The governing system is written in conservative form as

```math
\frac{\partial \mathbf{U}}{\partial t}
+
\nabla \cdot \mathbf{F}_c
-
\nabla \cdot \mathbf{F}_v
=
\mathbf{S}
```

where $\mathbf{F}_c$ denotes convective flux, $\mathbf{F}_v$ viscous and thermal flux, and $\mathbf{S}$ the axisymmetric geometric source contribution.

---

## Mass Conservation

```math
\frac{\partial \rho}{\partial t}
+
\nabla \cdot
\left(
\rho \mathbf{u}
\right)
=
0
```

---

## Momentum Conservation

```math
\frac{\partial (\rho \mathbf{u})}{\partial t}
+
\nabla \cdot
\left(
\rho \mathbf{u}\otimes\mathbf{u}
+
p\mathbf{I}
-
\boldsymbol{\tau}
\right)
=
\mathbf{S}_m
```

where $p$ is static pressure and $\boldsymbol{\tau}$ is the Newtonian viscous stress tensor.

---

## Energy Conservation

```math
\frac{\partial (\rho E)}{\partial t}
+
\nabla \cdot
\left[
(\rho E+p)\mathbf{u}
-
\boldsymbol{\tau}\cdot\mathbf{u}
+
\mathbf{q}
\right]
=
S_E
```

Thermal conduction follows Fourier's law:

```math
\mathbf{q}
=
-k\nabla T
```

---

## Thermodynamic Closure

V1 assumes a calorically perfect ideal gas:

```math
p
=
(\gamma-1)\rho e
```

The local speed of sound is

```math
a
=
\sqrt{\frac{\gamma p}{\rho}}
```

and Mach number is

```math
M
=
\frac{\lVert\mathbf{u}\rVert}{a}
```

Thermal conductivity is related to viscosity through the Prandtl number:

```math
k
=
\frac{\mu c_p}{Pr}
```

The current model deliberately excludes reacting chemistry and detailed multi-species rocket-exhaust thermodynamics.

---

# Numerical Formulation

## Conservative Finite Volumes

AstraFlow uses a **cell-centred finite-volume discretisation**.

For control volume $V_i$:

```math
\frac{d\mathbf{U}_i}{dt}
=
-
\frac{1}{V_i}
\sum_f
\mathbf{F}_f A_f
+
\mathbf{S}_i
```

The formulation evolves conserved quantities directly and is suitable for compressible flows containing shocks, contact discontinuities, rarefactions, and transonic acceleration.

---

## MUSCL Reconstruction

Second-order spatial reconstruction is performed using **MUSCL**.

Implemented slope limiters include:

- Minmod;
- Van Leer;
- Monotonized Central (MC).

The limiter controls nonphysical oscillations around steep gradients while preserving higher-order behaviour in smooth regions.

---

## HLLC / HLLE Fluxes

Convective intercell fluxes use the **HLLC approximate Riemann solver**.

HLLC resolves the approximate left, contact, and right wave families relevant to compressible Euler flow.

A more diffusive **HLLE fallback** is used when an HLLC reconstruction would otherwise produce a pathological intermediate state.

Fallback usage is exposed diagnostically rather than treated as an invisible normal operating mode.

---

## Time Integration

The production baseline uses **second-order Strong Stability Preserving Runge–Kutta integration (SSP-RK2)**.

The convective timestep is controlled by a CFL condition of the form

```math
\Delta t
\propto
\mathrm{CFL}
\frac{\Delta x}
{\lVert\mathbf{u}\rVert+a}
```

Viscous simulations additionally account for the appropriate diffusive stability restriction.

---

# Axisymmetric Geometry

AstraFlow models the meridional $(x,r)$ plane of a rotationally symmetric nozzle.

```text
r
↑
│
│  chamber
│  ──────────────────╲
│                     ╲
│                      ╲
│                       ╲ throat
│                        ╲
│                         ╲
│                          ╲──────────── exit
│
└──────────────────────────────────────→ x
```

The nozzle is parameterised by quantities including:

- chamber radius;
- chamber length;
- throat radius;
- contraction length;
- exit radius;
- expansion length;
- axial resolution;
- radial resolution.

A structured body-conforming grid follows the nozzle contour.

The finite-volume formulation incorporates axisymmetric geometry explicitly, including centreline treatment as $r \rightarrow 0$.

---

# Boundary Conditions

### Centreline

Axisymmetric symmetry at $r=0$.

### Inviscid Wall

Slip and impermeability for Euler calculations.

### Viscous Wall

No-slip, impermeable, adiabatic treatment for Navier–Stokes calculations.

### Inlet

Reservoir conditions are specified using stagnation quantities

```math
P_0,\qquad T_0
```

### Outlet

Back pressure is applied where the outlet remains subsonic.

For supersonic outflow, downstream quantities are not incorrectly imposed on characteristics that cannot propagate upstream.

---

# CPU and CUDA Backends

The CPU implementation acts as a numerical reference backend for:

- deterministic verification;
- debugging;
- regression testing;
- CPU-only CI;
- CPU/CUDA parity studies.

The CUDA implementation is the primary high-performance backend.

Development hardware:

| Property | Value |
|---|---|
| GPU | NVIDIA GeForce RTX 5070 Laptop GPU |
| Architecture | NVIDIA Blackwell |
| Compute capability | `12.0` |
| CUDA target | `sm_120` |

A project-local CUDA 12.8 toolchain has successfully compiled and executed native `sm_120` kernels on the target GPU.

---

# CUDA Architecture

Primary simulation fields use a **Structure of Arrays (SoA)** representation:

```text
rho[]
rho_u[]
rho_v[]
rho_E[]

pressure[]
temperature[]
mach[]

gradients[]
fluxes[]
residuals[]
```

The solver pipeline is conceptually:

```text
primitive variables
        ↓
gradients
        ↓
MUSCL reconstruction
        ↓
HLLC / HLLE fluxes
        ↓
viscous fluxes
        ↓
axisymmetric source terms
        ↓
residual assembly
        ↓
SSP-RK update
        ↓
residual reduction
```

Reusable device buffers avoid repeated allocation inside the timestep loop.

---

# Precision

AstraFlow supports FP32 and FP64 execution where required.

**FP32** is the primary interactive/high-throughput CUDA mode.

**FP64** is used extensively for numerical verification and CPU/GPU parity analysis.

Dimensional problems are internally scaled where appropriate to improve floating-point conditioning.

---

# Verification Results

The following values are **measured results from completed Phase 1 milestone tests**.

## Sod Shock Tube

For a 400-cell Sod problem:

| Quantity | Mean absolute error |
|---|---:|
| Density | `1.446 × 10^-3` |
| Velocity | `2.546 × 10^-3` |
| Pressure | `8.91 × 10^-4` |

Measured mass-conservation error:

```math
3.33 \times 10^{-15}
```

<!--
<p align="center">
  <img src="docs/assets/sod-verification.png" width="90%" alt="Sod verification">
</p>
-->

---

## CPU / CUDA Euler Parity

| Precision | Maximum CPU/GPU difference |
|---|---:|
| FP64 | `4.22 × 10^-14` |
| FP32 | `4.62 × 10^-6` |

Both precision paths passed their configured numerical parity criteria.

---

## 2D Conservative Transport

All four integrated conservative quantities were preserved within approximately

```math
10^{-12}
```

Maximum CPU/GPU difference:

```math
1.78 \times 10^{-15}
```

---

## Axisymmetric Nozzle Verification

| Metric | Result |
|---|---:|
| Throat Mach number | `1.0096` |
| Mean Mach difference from quasi-1D relation | `3.19 × 10^-4` |
| Axial mass-flow spread | `0.0774 %` |
| Stationary-state maximum drift | `2.22 × 10^-15` |
| CPU/GPU nozzle difference | `4.88 × 10^-15` |

The computed throat state is consistent with the expected transition through approximately sonic flow.

<!--
<p align="center">
  <img src="docs/assets/nozzle-verification.png" width="90%" alt="Nozzle verification">
</p>
-->

---

## Viscous Verification

Compressible Couette-flow verification currently gives:

| Quantity | Mean error |
|---|---:|
| Velocity | `5.15 × 10^-5` |
| Temperature | `4.46 × 10^-5` |

CPU/GPU parity:

| Precision | Maximum difference |
|---|---:|
| FP64 | `4.00 × 10^-15` |
| FP32 | `4.26 × 10^-6` |

---

## Final V1 Verification

| Check | Result |
|---|---:|
| Smooth entropy-wave refinement, CPU FP64 | L1 orders 1.81 and 2.00 |
| Couette temperature refinement, CPU FP64 | L1 orders 2.00 and 2.00 |
| Nozzle recovery from 20% slower initial flow | Throat Mach 0.803 → 0.997 |
| Long-time Couette CUDA FP32 | Velocity L1 5.94e-5; temperature L1 3.44e-5 |
| Release CPU / CUDA CTest | 6/6 and 8/8 passed |
| CPU ASan / UBSan / leak detection | 6/6 passed |
| WSLg GUI controls and real rendering | Passed scripted smoke/capture |

Both 128x32 rocket CLI runs completed 3,000 steps with valid JSON/CSV/VTK output and zero flux/reconstruction fallbacks. CUDA FP32 estimated mass flow 0.09661684 kg/s, thrust 171.05877 N and inlet/exit mass-flow mismatch 0.04655%. These runs stopped at the iteration limit, not the 1e-6 residual target; their outputs are transient ideal-gas estimates.

See the [verification report](docs/verification/verification_report.md) for analytical definitions, preset bounds, precision comparisons, grid tables and the reproducible sanitizer command. Grid refinement verifies the canonical cases; it does not establish grid independence of the coarse rocket example. Compute-sanitizer and Nsight were unavailable, so no passes from those tools are claimed.

---

# Engineering Analysis

## Mass Flow

```math
\dot{m}
=
\int_A \rho u\,dA
```

Mass-flow consistency is evaluated at multiple axial stations as a conservation diagnostic.

---

## Exit State

The numerical solution is integrated to obtain representative:

- Mach number;
- axial velocity;
- pressure;
- temperature.

Area- or mass-weighted averages are used where appropriate.

---

## Thrust

Nozzle thrust is obtained from the computed exit plane:

```math
F
=
\int_{A_e}
\rho u^2\,dA
+
\int_{A_e}
(p-p_a)\,dA
```

This permits nonuniform exit profiles to contribute directly to the result.

---

## Specific Impulse

```math
I_{sp}
=
\frac{F}{\dot{m}g_0}
```

These quantities describe the idealised numerical model rather than certification-grade engine performance.

---

# Reproducible Output

Each production simulation can generate:

```text
runs/<run-id>/
├── config.json
├── convergence.csv
├── performance.json
├── summary.json
├── mesh.vts
└── final_state.vts
```

| File | Purpose |
|---|---|
| `config.json` | effective simulation configuration |
| `convergence.csv` | residual history |
| `performance.json` | runtime/performance measurements |
| `summary.json` | engineering quantities and metadata |
| `mesh.vts` | structured computational mesh |
| `final_state.vts` | final CFD state |

VTK StructuredGrid output permits independent inspection in software such as ParaView. Explicit output paths must be new or empty. `scripts/check_run.py` validates the JSON/CSV/VTK output using only Python's standard library. Residuals are RMS time derivatives in nondimensional solver units; fields, times and engineering results use physical configuration units (SI in the rocket example). Each run explicitly distinguishes convergence, end-time and iteration-limit termination.

---

# Interfaces

## CLI

Example:

```bash
./build/astraflow_cli \
    --config examples/rocket_nozzle/config.json \
    --backend cuda
```

`--output runs/my-rocket`, `--max-iterations 3000` and `--precision float|double` override configuration values. Examples cover [Sod](examples/sod/config.json), an [isentropic nozzle](examples/isentropic_nozzle/config.json), a [viscous Couette channel](examples/viscous_channel/config.json) and a [viscous rocket nozzle](examples/rocket_nozzle/config.json).

Supported solver backends:

```text
cpu
cuda
```

The CLI reports quantities including:

- backend;
- mesh;
- iteration count;
- simulated time;
- wall-clock time;
- mass flow;
- throat Mach number;
- exit Mach number;
- exit velocity;
- exit pressure;
- thrust;
- specific impulse;
- mass-conservation error;
- residual.

---

## Interactive GUI

The V1 application uses:

- Dear ImGui;
- ImPlot;
- GLFW;
- OpenGL.

It invokes the same numerical solver as the CLI.

Controls include:

- Run;
- Pause;
- Single Step;
- Reset;
- Regenerate Mesh.

Visualisable quantities include:

- pressure;
- density;
- temperature;
- Mach number;
- axial velocity;
- radial velocity;
- velocity magnitude;
- total energy;
- vorticity.

Additional panels expose:

- convergence residuals;
- chamber conditions;
- nozzle geometry;
- engineering outputs;
- CUDA device information;
- iteration performance.

Simulation execution is isolated from interface rendering through a worker/state architecture. Immutable snapshots retain their own gas properties and remain valid across reset. Pause before applying geometry changes. Mouse-wheel zoom and middle-button pan change only the view.

Actual WSLg startup, OpenGL rendering, all nine fields, Run/Pause/Step/Reset/Regenerate and residual/engineering updates passed `--smoke-test --capture runs/gui.ppm`. This is scripted control-handler testing and inspected framebuffer output, not manual mouse testing.

---

# Performance Study

Release measurements compare equal FP32 on the RTX 5070 Laptop GPU and a single-thread CPU reference on an Intel Core Ultra 9 285H. Each grid uses 20 warmup iterations and 50 measured iterations; setup and output are excluded.

| Grid | CPU ms/iteration | GPU ms/iteration | CPU/GPU | Device MiB |
|---|---:|---:|---:|---:|
| 128 × 32 | 1.448 | 0.539 | 2.69x | 1.54 |
| 256 × 64 | 6.194 | 0.512 | 12.11x | 6.10 |
| 512 × 128 | 28.222 | 0.619 | 45.56x | 24.32 |
| 1024 × 256 | 121.539 | 1.987 | 61.16x | 97.11 |

These are measured sample means, not guaranteed speedups or comparisons against a tuned multicore CPU solver. Small problems can be CPU-faster: the earlier 400-cell 1D case was slower on GPU. Laptop clock/load variation also affected the before/after measurements.

CUDA-event profiling identified diagnostic reductions/transfers as the largest measured stage. A single four-component CUB reduction removed three reduction calls and the residual-square buffer. The measured diagnostic interval at 512x128 fell from 0.306395 to 0.184032 ms, and numerical parity still passes. Instrumentation adds overhead; uninstrumented iteration times are reported above. Buffer counts exclude driver/context and OpenGL resources.

Full methodology, limitations, raw before/after measurements and reproduction commands are in the [benchmark report](docs/performance/benchmark_report.md). Enable the benchmark executable with `-DASTRAFLOW_BUILD_BENCHMARKS=ON`.

---

# Build

## Requirements

```text
CMake >= 3.24
Ninja
C++20-compatible compiler
```

CUDA builds require CUDA 12.8+ with native SM120 support. Tested with GCC 13.3, CMake 3.28.3 and CUDA 12.8.93 on Ubuntu 24.04 / WSL2. CMake downloads pinned dependencies at first configure. No global Python packages are needed.

Primary development target:

```text
NVIDIA Blackwell
compute_120
sm_120
```

## CPU

```bash
cmake \
    -S . \
    -B build-cpu \
    -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DASTRAFLOW_ENABLE_CUDA=OFF

cmake --build build-cpu
```

## CUDA

```bash
cmake \
    -S . \
    -B build \
    -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DASTRAFLOW_ENABLE_CUDA=ON

cmake --build build
```

## Project-local CUDA and GUI

If CUDA is not already installed, the bootstrap downloads SHA256-checked compiler/runtime/CCCL packages into this project only:

```bash
python3 -m venv .venv
.venv/bin/python scripts/bootstrap_cuda.py
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CUDA_COMPILER="$PWD/.toolchains/cuda-12.8.1/bin/nvcc" \
    -DASTRAFLOW_ENABLE_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=120
cmake --build build -j 6
./build/astraflow_cli --device-info
```

The device probe executes a kernel and reports its compiled architecture. Agent sandbox restrictions may require GPU execution outside the sandbox.

```bash
# Ubuntu 24.04 only, if GUI development headers are missing:
bash scripts/bootstrap_gui.sh
cmake -S . -B build -DASTRAFLOW_BUILD_GUI=ON
cmake --build build -j 6
./build/astraflow_gui --config examples/rocket_nozzle/config.json
```

The GUI bootstrap extracts development headers/libraries project-locally and reuses the system OpenGL/X11 runtime. `ASTRAFLOW_BUILD_TESTS` defaults ON; GUI and benchmarks default OFF. CPU-only GitHub Actions are configured; GPU verification runs locally on the target hardware.

## Tests

```bash
ctest \
    --test-dir build \
    --output-on-failure
```

CPU-only:

```bash
ctest \
    --test-dir build-cpu \
    --output-on-failure
```

---

# Technology

| Area | Technology |
|---|---|
| Core implementation | C++20 |
| GPU computing | CUDA C++ |
| Build | CMake + Ninja |
| Configuration | nlohmann/json |
| Verification | Catch2 |
| GUI | Dear ImGui |
| Plotting | ImPlot |
| Windowing | GLFW |
| Rendering | OpenGL |
| Scientific output | JSON / CSV / VTK StructuredGrid |
| CI | GitHub Actions |

Pinned dependencies:

- `nlohmann/json 3.11.3` (MIT);
- `Catch2 3.7.1` (BSL-1.0);
- `Dear ImGui 1.91.8` (MIT);
- `ImPlot 0.16` (MIT);
- `GLFW 3.4` (zlib/libpng).

Their upstream license files remain in fetched sources. Auxiliary Python scripts use only the standard library.

---

# Architecture

```text
AstraFlow/
├── apps/
│   ├── cli/
│   └── gui/
├── include/
│   └── astraflow/
├── src/
│   ├── core/
│   ├── cpu/
│   ├── cuda/
│   ├── physics/
│   ├── numerics/
│   ├── geometry/
│   ├── analysis/
│   ├── io/
│   └── visualization/
├── tests/
│   ├── unit/
│   ├── verification/
│   └── regression/
├── benchmarks/
├── examples/
├── scripts/
├── docs/
└── .github/
```

Dependency direction is intentionally constrained:

```text
Mathematics / physics
        ↓
Numerical methods
        ↓
CPU / CUDA backends
        ↓
Analysis / I/O
        ↓
CLI / GUI
```

**Numerical physics is not implemented inside GUI code.**

See the [architecture overview](docs/architecture/overview.md), [governing equations](docs/mathematics/governing_equations.md), [numerical method](docs/mathematics/numerical_method.md), [axisymmetric formulation](docs/mathematics/axisymmetric_navier_stokes.md), [boundary conditions](docs/mathematics/boundary_conditions.md) and [environment record](docs/environment.md).

---

# Development Milestones

| Milestone | Scope | Status |
|---|---|---|
| M0 | Toolchain, repository, CUDA `sm_120` | Complete |
| M1 | 1D CPU Euler reference solver | Complete |
| M2 | CUDA Euler backend | Complete |
| M3 | 2D finite-volume infrastructure | Complete |
| M4 | Axisymmetric nozzle solver | Complete |
| M5 | Viscous compressible Navier–Stokes | Complete |
| M6 | Engineering analysis, CLI and output | Complete |
| M7 | Interactive GUI | Complete |
| M8 | CUDA profiling and optimisation | Complete |
| M9 | Final implementation verification | Complete |

---

# V1 Scope

V1 deliberately excludes:

- reacting combustion chemistry;
- multi-species transport;
- turbulence modelling;
- LES;
- full 3D geometry;
- external exhaust plumes;
- adaptive mesh refinement;
- conjugate wall heat transfer;
- regenerative cooling;
- real-gas thermodynamics;
- multi-GPU execution.

These remain candidates for later development after V1 verification.

---

# Planned Extensions

Potential post-V1 directions include:

- external plume simulation;
- temperature-dependent thermodynamics;
- multi-species reacting flow;
- turbulence modelling;
- wall heat-flux analysis;
- conjugate heat transfer;
- regenerative cooling;
- CUDA/OpenGL interoperability;
- additional GPU optimisation;
- nozzle-geometry optimisation;
- larger and eventually three-dimensional simulations.

Future functionality will be introduced alongside appropriate verification cases.

---

# Scientific Integrity

Project-level requirements include:

- quantitative verification rather than visual validation;
- CPU/CUDA cross-checking;
- physically justified tolerances;
- explicit invalid-state detection;
- measured rather than assumed benchmark results;
- documented physical limitations;
- preservation of verification after optimisation;
- no claim of fidelity beyond the implemented physical model.

AstraFlow numerically solves a discretised compressible Euler/Navier–Stokes model. It makes no claim regarding the mathematical Navier–Stokes existence-and-smoothness problem.

---

# Author

**Arnav**  
BSc Computer Science with Artificial Intelligence  
University of Nottingham

Project responsibilities:

- mathematical formulation;
- numerical-method architecture;
- CFD system design;
- CUDA architecture;
- verification strategy;
- software architecture;
- technical review;
- direction of AI-assisted implementation.

---

# License

AstraFlow is released under the **MIT License**.

Third-party libraries retain their respective upstream licences. NVIDIA CUDA and associated tooling remain subject to NVIDIA's applicable licence terms.

---

# Disclaimer

AstraFlow is an educational, research, and portfolio scientific-computing project.

It is not flight-qualified aerospace software and should not be used as the sole basis for safety-critical aerospace design, manufacture, or operation.

Results must be interpreted within the assumptions, model fidelity, numerical discretisation, and verification evidence documented by the project.

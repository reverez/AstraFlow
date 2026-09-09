# AstraFlow

> **GPU-accelerated compressible CFD for axisymmetric rocket nozzle flows.**

AstraFlow is a scientific-computing project implementing a **2D axisymmetric compressible Euler/Navier–Stokes solver** in **C++20 and NVIDIA CUDA**, with a CPU reference backend, scientific verification suite, engineering analysis pipeline, reproducible simulation output, and interactive visualisation.

The project targets an **NVIDIA GeForce RTX 5070 Laptop GPU** using NVIDIA Blackwell **compute capability 12.0 (`sm_120`)**.

I architected the mathematical formulation, numerical method, physical assumptions, verification strategy, CPU/GPU architecture, project structure, and milestone acceptance criteria. AI coding agents are used to accelerate implementation, testing, debugging, profiling, and documentation under those specifications.

Numerical results are accepted only after comparison against analytical solutions, conservation properties, convergence behaviour, or CPU/GPU parity tests.

---

## Status

**Phase 1 / V1 is under active development.**

Milestones **0–6 are complete**, including:

- verified 1D CPU Euler solver;
- CUDA Euler backend;
- 2D conservative finite-volume infrastructure;
- axisymmetric rocket-nozzle formulation;
- viscous compressible Navier–Stokes transport;
- CPU/CUDA numerical parity;
- engineering analysis;
- headless CLI simulation;
- JSON, CSV, and VTK output.

Remaining work consists primarily of:

- interactive GUI completion;
- CUDA profiling and optimisation;
- grid-refinement verification;
- final CPU/GPU benchmarking;
- V1 acceptance review.

See [`docs/PHASE1_STATUS.md`](docs/PHASE1_STATUS.md).

---

## Preview

<!-- Replace these placeholders as Phase 1 is completed. -->

### Interactive Simulation

<!--
<p align="center">
  <img src="docs/assets/gui-overview.png" width="100%" alt="AstraFlow interactive CFD interface">
</p>
-->

> **Image pending:** interactive nozzle simulation, convergence history, engineering analysis, and CUDA telemetry.

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

## Remaining V1 Verification

Before V1 acceptance:

- grid-refinement analysis;
- final Release CPU regression;
- final Release CUDA regression;
- full rocket-nozzle CPU/CUDA comparison;
- CUDA memory/error diagnostics;
- GUI acceptance;
- final performance profiling.

Formal convergence claims will not be made until the grid-refinement study is complete.

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

VTK StructuredGrid output permits independent inspection in software such as ParaView.

---

# Interfaces

## CLI

Example:

```bash
astraflow_cli \
    --config examples/rocket_nozzle/config.json \
    --backend cuda
```

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

Simulation execution is isolated from interface rendering through a worker/state architecture.

---

# Performance Study

The final V1 profiling pass will evaluate representative mesh sizes such as:

```text
128 × 32
256 × 64
512 × 128
1024 × 256
```

subject to practical memory and runtime constraints.

Metrics include:

| Metric | Purpose |
|---|---|
| Cell count | problem scale |
| CPU iteration time | reference performance |
| CUDA iteration time | GPU performance |
| Iterations/s | solver throughput |
| CPU/GPU speedup | acceleration |
| GPU memory | memory scaling |
| Kernel timings | profiling |

No minimum CUDA speedup is assumed in advance.

Small meshes may remain CPU-faster because GPU launch and synchronisation overhead can dominate limited parallel work. Final claims will use measured scaling data only.

<!--
<p align="center">
  <img src="docs/assets/cpu-gpu-scaling.png" width="90%" alt="CPU versus CUDA scaling">
</p>
-->

---

# Build

## Requirements

```text
CMake >= 3.24
Ninja
C++20-compatible compiler
```

CUDA builds additionally require an NVIDIA CUDA toolchain compatible with the configured device architecture.

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

Current pinned core dependencies include:

- `nlohmann/json 3.11.3`;
- `Catch2 3.7.1`.

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
| M7 | Interactive GUI | In progress |
| M8 | CUDA profiling and optimisation | Pending |
| M9 | Final V1 acceptance | Pending |

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

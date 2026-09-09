# AstraFlow

> **GPU-accelerated compressible CFD for axisymmetric rocket nozzle flows.**

AstraFlow is a scientific-computing project implementing a **2D axisymmetric compressible Euler/Navier–Stokes solver** in **C++20 and NVIDIA CUDA**, with a CPU reference backend, scientific verification suite, engineering analysis pipeline, reproducible simulation output, and interactive visualisation.

The project is architected around a conservative finite-volume formulation for high-speed compressible flow and targets an **NVIDIA GeForce RTX 5070 Laptop GPU** using native NVIDIA Blackwell **compute capability 12.0 (`sm_120`)**.

I designed the mathematical model, numerical architecture, physical assumptions, verification strategy, GPU/CPU division, project structure, and development milestones. AI coding agents are used to accelerate implementation, testing, debugging, and documentation under those specifications. Numerical results are accepted only after comparison against analytical solutions, conservation properties, or CPU/GPU parity tests.

---

## Status

**Phase 1 / V1 is under active development.**

Milestones **0–6 are complete**, covering:

- verified 1D CPU Euler solver;
- CUDA Euler solver;
- 2D finite-volume infrastructure;
- axisymmetric nozzle formulation;
- viscous compressible Navier–Stokes transport;
- CPU/CUDA numerical parity;
- engineering analysis;
- command-line simulation;
- JSON, CSV, and VTK output.

The interactive GUI, CUDA optimisation pass, grid-refinement study, final benchmark suite, and V1 acceptance review remain in progress.

See [`docs/PHASE1_STATUS.md`](docs/PHASE1_STATUS.md) for the implementation ledger.

---

## Preview

<!-- Replace these placeholders as Phase 1 visualisation and benchmarking are completed. -->

### Interactive Simulation

<!--
<p align="center">
  <img src="docs/assets/gui-overview.png" width="100%" alt="AstraFlow interactive CFD interface">
</p>
-->

> **Image pending:** AstraFlow GUI showing the nozzle flow field, simulation controls, residual history, engineering quantities, and CUDA performance.

### Mach Field

<!--
<p align="center">
  <img src="docs/assets/nozzle-mach.png" width="95%" alt="AstraFlow Mach-number field">
</p>
-->

> **Image pending:** Mach-number contours through the converging-diverging nozzle.

### Pressure Field

<!--
<p align="center">
  <img src="docs/assets/nozzle-pressure.png" width="95%" alt="AstraFlow pressure field">
</p>
-->

> **Image pending:** Static-pressure distribution through the nozzle.

### CPU vs CUDA Scaling

<!--
<p align="center">
  <img src="docs/assets/cpu-gpu-scaling.png" width="90%" alt="AstraFlow CPU versus CUDA performance">
</p>
-->

> **Graph pending:** measured CPU/CUDA iteration time and speedup across increasing mesh sizes.

---

# Project Objectives

AstraFlow is intended to investigate the complete path from continuous mathematical model to verified high-performance numerical software.

The core objectives are to:

1. implement the governing equations rather than rely on an external CFD solver;
2. construct a conservative numerical method suitable for compressible and transonic flow;
3. reproduce characteristic internal rocket-nozzle behaviour;
4. verify the numerical implementation against known solutions and invariants;
5. maintain independent CPU and CUDA execution paths for cross-validation;
6. exploit a modern consumer NVIDIA GPU for parallel finite-volume computation;
7. quantify numerical error and hardware performance rather than relying on visual plausibility;
8. expose the solver through both a reproducible CLI and an interactive engineering interface.

The project is therefore both a **CFD solver** and a study in **numerical methods, GPU computing, scientific software engineering, and AI-assisted technical development**.

---

# Authorship and AI-Assisted Development

AstraFlow is developed using an **architect-directed, AI-assisted implementation workflow**.

I am responsible for defining and reviewing:

- physical scope;
- governing equations;
- state representation;
- finite-volume formulation;
- reconstruction method;
- Riemann solver;
- temporal integration;
- axisymmetric treatment;
- viscous and thermal transport;
- boundary conditions;
- numerical stability requirements;
- verification cases;
- error tolerances;
- CPU/CUDA comparison strategy;
- GPU architecture;
- software architecture;
- milestone acceptance criteria;
- interpretation of numerical results.

AI coding agents, principally Codex, are used as implementation tools to accelerate:

- C++ and CUDA implementation;
- test construction;
- build-system work;
- debugging;
- numerical experiments;
- repetitive infrastructure;
- profiling;
- documentation.

The development process is not based on accepting generated code because it compiles or produces plausible images.

Each numerical subsystem is required to satisfy targeted tests before becoming part of the accepted solver:

```text
Mathematical design
        ↓
Implementation specification
        ↓
AI-assisted implementation
        ↓
Targeted numerical tests
        ↓
Analytical / conservation / parity validation
        ↓
Architectural review
        ↓
Milestone acceptance
```

This separation is particularly important for scientific computing, where an implementation can execute successfully while remaining physically or numerically incorrect.

---

# Physical Model

V1 models a **single-species, calorically perfect ideal gas** flowing through an axisymmetric converging-diverging nozzle.

The conservative state is

$$
\mathbf{U}
=
\begin{bmatrix}
\rho \\
\rho u \\
\rho v \\
\rho E
\end{bmatrix},
$$

where:

- $\rho$ is density;
- $u$ is axial velocity;
- $v$ is radial velocity;
- $E$ is total specific energy.

The governing system can be written in conservative form as

$$
\frac{\partial \mathbf{U}}{\partial t}
+
\nabla \cdot \mathbf{F}_c
-
\nabla \cdot \mathbf{F}_v
=
\mathbf{S},
$$

where $\mathbf{F}_c$ contains convective fluxes, $\mathbf{F}_v$ contains viscous and thermal fluxes, and $\mathbf{S}$ represents axisymmetric geometric source terms.

---

## Mass Conservation

$$
\frac{\partial \rho}{\partial t}
+
\nabla \cdot (\rho \mathbf{u})
=
0.
$$

---

## Momentum Conservation

$$
\frac{\partial (\rho \mathbf{u})}{\partial t}
+
\nabla \cdot
\left(
\rho \mathbf{u} \otimes \mathbf{u}
+
p\mathbf{I}
-
\boldsymbol{\tau}
\right)
=
\mathbf{S}_m.
$$

Here $p$ is static pressure and $\boldsymbol{\tau}$ is the Newtonian viscous stress tensor.

---

## Energy Conservation

$$
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
S_E.
$$

Thermal conduction follows Fourier's law,

$$
\mathbf{q}
=
-k\nabla T.
$$

---

## Thermodynamics

V1 uses the calorically perfect ideal-gas relation

$$
p
=
(\gamma-1)\rho e,
$$

with local speed of sound

$$
a
=
\sqrt{\frac{\gamma p}{\rho}},
$$

and Mach number

$$
M
=
\frac{\lVert \mathbf{u} \rVert}{a}.
$$

Thermal conductivity is related to viscosity using the Prandtl number,

$$
k
=
\frac{\mu c_p}{Pr}.
$$

The model intentionally does not represent detailed rocket combustion chemistry or chemically reacting exhaust.

---

# Numerical Formulation

## Finite-Volume Discretisation

AstraFlow uses a **cell-centred conservative finite-volume scheme**.

For control volume $V_i$,

$$
\frac{d\mathbf{U}_i}{dt}
=
-
\frac{1}{V_i}
\sum_f
\mathbf{F}_f A_f
+
\mathbf{S}_i.
$$

The formulation directly evolves conserved quantities and is suitable for compressible flow containing strong gradients, shocks, contact discontinuities, and transonic acceleration.

---

## Spatial Reconstruction

Second-order spatial reconstruction is performed using **MUSCL**.

Implemented slope limiters include:

- Minmod;
- Van Leer;
- Monotonized Central (MC).

These limit reconstructed gradients near discontinuities to suppress nonphysical oscillation while retaining higher-order behaviour in smooth regions.

---

## Convective Flux

Intercell convective fluxes use the **HLLC approximate Riemann solver**.

HLLC resolves the principal left-, contact-, and right-moving wave structure required for compressible Euler flow.

A more diffusive **HLLE fallback** is available when reconstructed states would otherwise produce a pathological intermediate solution.

Fallback usage is treated as a diagnostic rather than silently becoming the normal numerical method.

---

## Time Integration

The production baseline uses **second-order Strong Stability Preserving Runge–Kutta integration (SSP-RK2)**.

The timestep is controlled by a CFL stability condition based on local convective signal speed,

$$
\Delta t
\propto
\mathrm{CFL}
\frac{\Delta x}
{\lVert \mathbf{u} \rVert + a}.
$$

Viscous simulations additionally account for the appropriate diffusive stability restriction.

---

# Axisymmetric Formulation

AstraFlow models the meridional $(x,r)$ plane of a rotationally symmetric nozzle.

The computational domain therefore represents a three-dimensional axisymmetric geometry without requiring a full 3D mesh.

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

The nozzle geometry is parameterised by quantities including:

- chamber radius;
- chamber length;
- throat radius;
- contraction length;
- exit radius;
- expansion length;
- axial resolution;
- radial resolution.

A structured body-conforming mesh is generated from the nozzle contour.

Axisymmetric finite-volume geometry and source terms are included explicitly. The centreline treatment is constructed to avoid singular numerical evaluation as $r \rightarrow 0$.

---

# Boundary Conditions

V1 includes the boundary conditions required for internal nozzle simulation.

### Axis

Axisymmetric symmetry at $r=0$.

### Inviscid Wall

Slip and impermeability for Euler verification.

### Viscous Wall

No-slip, impermeable, and adiabatic treatment for Navier–Stokes simulations.

### Inlet

Reservoir conditions are specified using stagnation pressure and temperature,

$$
P_0,\qquad T_0.
$$

### Outlet

Back pressure is applied where the outlet remains subsonic.

Supersonic outflow is treated without incorrectly imposing all downstream primitive variables on a region from which information cannot propagate upstream.

---

# CPU and CUDA Backends

AstraFlow contains separate CPU and CUDA execution backends operating on the same mathematical model.

The CPU implementation serves primarily as a **reference backend** for:

- deterministic verification;
- debugging;
- regression testing;
- CPU-only CI;
- CPU/CUDA parity analysis.

The CUDA implementation is the production high-performance backend.

Primary development hardware:

```text
GPU:                  NVIDIA GeForce RTX 5070 Laptop GPU
Architecture:         NVIDIA Blackwell
Compute capability:   12.0
CUDA target:          sm_120
```

The project-local CUDA toolchain successfully compiles and executes native `sm_120` kernels on the target GPU.

---

# CUDA Architecture

Simulation fields use a **Structure of Arrays (SoA)** layout where appropriate:

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

This layout is intended to improve memory coalescing when neighbouring CUDA threads access the same physical field.

The GPU backend is structured around stages conceptually equivalent to:

```text
primitive-variable evaluation
        ↓
gradient reconstruction
        ↓
MUSCL face reconstruction
        ↓
HLLC / HLLE flux evaluation
        ↓
viscous flux evaluation
        ↓
axisymmetric source terms
        ↓
residual assembly
        ↓
SSP-RK update
        ↓
residual reduction
```

Device allocations are managed through reusable owned buffers rather than repeatedly allocating memory inside the timestep loop.

---

# Floating-Point Precision

The solver supports FP32 and FP64 execution where required.

### FP32

Primary interactive and high-throughput CUDA mode.

### FP64

Used extensively for verification and CPU/GPU parity analysis.

Physical problems are internally scaled where appropriate to improve conditioning and avoid unnecessary loss of precision when using dimensional rocket-engine quantities.

---

# Verification Results

The following values were measured during completed Phase 1 milestone tests.

They are not projected values.

## Sod Shock Tube

A 400-cell Sod problem produced:

| Quantity | Mean absolute error |
|---|---:|
| Density | `1.446 × 10^-3` |
| Velocity | `2.546 × 10^-3` |
| Pressure | `8.91 × 10^-4` |

Measured mass-conservation error:

$$
3.33 \times 10^{-15}.
$$

<!-- Add Sod comparison graph here. -->

<!--
<p align="center">
  <img src="docs/assets/sod-verification.png" width="90%" alt="Sod shock-tube numerical versus exact solution">
</p>
-->

---

## CPU / CUDA Euler Parity

Maximum measured CPU/GPU state difference:

| Precision | Maximum difference |
|---|---:|
| FP64 | `4.22 × 10^-14` |
| FP32 | `4.62 × 10^-6` |

Both passed their configured numerical parity criteria.

---

## 2D Conservative Transport

The transported-density-wave verification conserved all four integrated conservative quantities within approximately

$$
10^{-12}.
$$

Maximum measured CPU/GPU difference:

$$
1.78 \times 10^{-15}.
$$

---

## Axisymmetric Nozzle

The current inviscid nozzle verification produced:

| Metric | Result |
|---|---:|
| Throat Mach number | `1.0096` |
| Mean Mach difference from quasi-1D relation | `3.19 × 10^-4` |
| Axial mass-flow spread | `0.0774 %` |
| Stationary-state maximum drift | `2.22 × 10^-15` |
| CPU/GPU nozzle difference | `4.88 × 10^-15` |

The throat result demonstrates the expected transition through approximately sonic conditions.

<!-- Add analytical-vs-numerical nozzle graph here. -->

<!--
<p align="center">
  <img src="docs/assets/nozzle-verification.png" width="90%" alt="Numerical and quasi-1D nozzle Mach comparison">
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

<!-- Add Couette analytical comparison here. -->

---

## Remaining V1 Verification

Before V1 acceptance, the project will additionally complete:

- multi-resolution grid-refinement analysis;
- final Release CPU regression;
- final Release CUDA regression;
- full rocket-nozzle CPU/CUDA comparison;
- CUDA memory/error diagnostics;
- final GUI acceptance;
- final performance profiling.

Formal convergence claims will not be made until the grid-refinement study is complete.

---

# Engineering Analysis

The CFD solution is converted into nozzle-performance quantities directly from the numerical field.

## Mass Flow

$$
\dot{m}
=
\int_A \rho u\,dA.
$$

Mass-flow consistency is also evaluated across multiple axial stations as a conservation diagnostic.

---

## Exit Conditions

The solver calculates representative exit quantities including:

- Mach number;
- axial velocity;
- pressure;
- temperature.

Area- or mass-weighted averaging is used where physically appropriate.

---

## Thrust

Nozzle thrust is evaluated from the computed exit plane using

$$
F
=
\int_{A_e} \rho u^2\,dA
+
\int_{A_e}(p-p_a)\,dA.
$$

This avoids assuming that the exit state is perfectly uniform.

---

## Specific Impulse

$$
I_{sp}
=
\frac{F}{\dot{m}g_0}.
$$

These quantities describe the idealised numerical model and are not presented as experimental or certification-grade rocket-engine predictions.

---

# Reproducible Output

Each production simulation can generate a self-contained run directory:

```text
runs/<run-id>/
├── config.json
├── convergence.csv
├── performance.json
├── summary.json
├── mesh.vts
└── final_state.vts
```

### `config.json`

Effective simulation configuration.

### `convergence.csv`

Residual history.

### `performance.json`

Execution and performance data.

### `summary.json`

Derived engineering quantities and run metadata.

### `mesh.vts`

Structured computational mesh.

### `final_state.vts`

Final CFD state for independent visualisation.

VTK StructuredGrid output allows results to be inspected independently in applications such as ParaView.

---

# Interfaces

## CLI

AstraFlow supports headless simulation:

```bash
astraflow_cli \
    --config examples/rocket_nozzle/config.json \
    --backend cuda
```

Supported backends include:

```text
cpu
cuda
```

The CLI reports quantities including:

```text
Backend
Grid
Iterations
Simulated time
Wall-clock time
Mass flow
Throat Mach
Exit Mach
Exit velocity
Exit pressure
Thrust
Specific impulse
Mass-conservation error
Final residual
```

---

## Interactive GUI

The V1 interface uses:

- Dear ImGui;
- ImPlot;
- GLFW;
- OpenGL.

The GUI uses the same underlying solver as the CLI rather than implementing separate physics.

Planned/implemented controls include:

- Run;
- Pause;
- Single Step;
- Reset;
- Regenerate Mesh.

Visualisable fields include:

- pressure;
- density;
- temperature;
- Mach number;
- axial velocity;
- radial velocity;
- velocity magnitude;
- total energy;
- vorticity.

The interface also exposes:

- residual convergence;
- nozzle geometry parameters;
- chamber conditions;
- engineering quantities;
- CUDA device information;
- solver iteration performance.

Simulation execution is separated from GUI rendering through a worker/state model so long-running iterations do not block the interface.

---

# Performance Study

A dedicated profiling milestone will evaluate representative grids such as:

```text
128 × 32
256 × 64
512 × 128
1024 × 256
```

subject to practical memory and runtime constraints.

The final report will measure:

| Metric | Purpose |
|---|---|
| Cell count | problem scale |
| CPU iteration time | reference performance |
| CUDA iteration time | GPU performance |
| Iterations/s | numerical throughput |
| GPU speedup | acceleration |
| GPU memory | memory scaling |
| Kernel timings | optimisation analysis |

Small problems are expected to expose CUDA launch and synchronisation overhead, while larger meshes provide progressively more parallel work.

No minimum CUDA speedup is assumed in advance. Final performance claims will use measured results only.

<!-- Add final benchmark graph here. -->

---

# Build

## Requirements

CPU builds require:

```text
CMake >= 3.24
Ninja
C++20-compatible compiler
```

CUDA builds additionally require an NVIDIA CUDA toolchain capable of targeting the selected GPU architecture.

The primary development configuration targets:

```text
NVIDIA Blackwell
compute_120
sm_120
```

---

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

---

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

---

## Tests

```bash
ctest \
    --test-dir build \
    --output-on-failure
```

CPU-only testing:

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
| Build system | CMake + Ninja |
| Configuration | nlohmann/json |
| Testing | Catch2 |
| GUI | Dear ImGui |
| Plotting | ImPlot |
| Windowing | GLFW |
| Rendering | OpenGL |
| Scientific output | JSON / CSV / VTK StructuredGrid |
| CI | GitHub Actions |

Current pinned core dependencies include:

- `nlohmann/json 3.11.3`;
- `Catch2 3.7.1`.

Third-party components retain their respective upstream licences.

---

# Architecture

The repository is organised around explicit separation of numerical and application concerns:

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

The main architectural boundaries are:

```text
Mathematics / physics
        ↓
Numerical methods
        ↓
CPU / CUDA backends
        ↓
Analysis and I/O
        ↓
CLI / GUI
```

**Physics and numerical algorithms are not implemented inside GUI code.**

---

# Development Milestones

| Milestone | Scope | Status |
|---|---|---|
| M0 | Toolchain, CMake, CUDA `sm_120`, repository | Complete |
| M1 | Verified 1D CPU Euler solver | Complete |
| M2 | CUDA Euler backend and CPU/GPU parity | Complete |
| M3 | Verified 2D finite-volume infrastructure | Complete |
| M4 | Axisymmetric rocket-nozzle solver | Complete |
| M5 | Viscous compressible Navier–Stokes | Complete |
| M6 | Engineering analysis, CLI and output | Complete |
| M7 | Interactive GUI | In progress |
| M8 | CUDA profiling and optimisation | Pending |
| M9 | Final V1 scientific acceptance | Pending |

---

# V1 Scope Boundary

AstraFlow V1 deliberately excludes:

- reacting combustion chemistry;
- multi-species transport;
- turbulence models;
- LES;
- full 3D geometry;
- external exhaust-plume simulation;
- adaptive mesh refinement;
- conjugate wall heat transfer;
- regenerative cooling;
- real-gas thermodynamics;
- multi-GPU execution.

These are candidate extensions after the baseline solver has completed numerical and architectural acceptance.

---

# Planned Extensions

Following V1, potential development directions include:

- external under-expanded and over-expanded plume simulation;
- temperature-dependent thermodynamics;
- multi-species reacting flow;
- turbulence modelling;
- wall heat-flux analysis;
- conjugate heat transfer;
- regenerative-cooling studies;
- CUDA/OpenGL interoperability;
- further CUDA kernel optimisation;
- nozzle geometry optimisation;
- larger and eventually three-dimensional simulations.

Future additions will be introduced only where they can be independently verified.

---

# Scientific Integrity

AstraFlow follows several project-level rules:

- numerical results must be verified rather than inferred from appearance;
- CUDA results must be compared against a reference implementation;
- tolerances must not be loosened solely to make tests pass;
- NaN, non-finite, negative-density, and negative-pressure states are treated as failures;
- benchmark results must be measured rather than estimated;
- numerical limitations must be documented;
- the ideal-gas nozzle model must not be presented as complete rocket-engine physics;
- performance optimisation must preserve verification results.

AstraFlow numerically solves a discretised compressible Euler/Navier–Stokes model. It makes no claim regarding the mathematical Navier–Stokes existence-and-smoothness problem.

---

# Author

**Arnav**  
BSc Computer Science with Artificial Intelligence  
University of Nottingham

Project roles:

- mathematical formulation;
- numerical architecture;
- CFD system design;
- CUDA architecture;
- verification strategy;
- software architecture;
- technical review;
- AI-assisted implementation direction.

---

# License

AstraFlow is released under the **MIT License**.

Third-party libraries retain their respective upstream licences. NVIDIA CUDA and associated tooling remain subject to NVIDIA's applicable licence terms.

---

# Disclaimer

AstraFlow is an educational, research, and portfolio scientific-computing project.

It is not flight-qualified software and should not be used as the sole basis for safety-critical aerospace design, manufacture, or operation.

Results must be interpreted within the numerical assumptions, discretisation error, model fidelity, and verification evidence documented by the project.

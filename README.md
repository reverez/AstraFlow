# AstraFlow

> **GPU-accelerated compressible computational fluid dynamics for axisymmetric rocket nozzle flows.**

AstraFlow is a scientific-computing and high-performance computing project that implements a **2D axisymmetric compressible Euler/Navier–Stokes solver** for internal rocket-nozzle flow using **C++20 and NVIDIA CUDA**.

The project is designed around two complementary solver backends:

- a **CPU reference implementation** for numerical verification, deterministic testing, and debugging;
- a **CUDA production implementation** designed to exploit massively parallel NVIDIA GPU hardware.

AstraFlow combines numerical fluid dynamics, thermodynamics, partial differential equations, finite-volume methods, GPU programming, scientific verification, software engineering, and interactive visualisation into a single portfolio-scale system.

The current target platform is an **NVIDIA GeForce RTX 5070 Laptop GPU**, using NVIDIA Blackwell compute capability **12.0 (`sm_120`)**.

---

> [!NOTE]
> **Development status — Phase 1 / V1**
>
> Core solver development is substantially complete. Milestones 0–6 have passed, covering the CPU and CUDA Euler solvers, 2D finite-volume infrastructure, axisymmetric rocket-nozzle formulation, viscous Navier–Stokes transport, engineering analysis, CLI execution, and reproducible scientific output.
>
> Interactive GUI integration, CUDA performance profiling, grid-refinement verification, final documentation, and the comprehensive V1 acceptance pass remain in progress.
>
> Detailed implementation status is maintained in [`docs/PHASE1_STATUS.md`](docs/PHASE1_STATUS.md).

---

## Preview

<!-- Replace these placeholders with final project screenshots after GUI completion. -->

### Interactive CFD Workspace

> **IMAGE PLACEHOLDER — AstraFlow GUI**
>
> Recommended final image:
> - full AstraFlow application window;
> - nozzle Mach-number field;
> - simulation controls;
> - convergence graph;
> - engineering metrics;
> - RTX 5070 / CUDA performance panel.

<!--
<p align="center">
  <img src="docs/assets/gui-overview.png" width="100%" alt="AstraFlow interactive CFD interface">
</p>
-->

### Rocket Nozzle Flow Field

> **IMAGE PLACEHOLDER — Mach / Pressure Contours**

<!--
<p align="center">
  <img src="docs/assets/nozzle-mach.png" width="90%" alt="AstraFlow rocket nozzle Mach field">
</p>
-->

### CPU vs CUDA Performance

> **IMAGE PLACEHOLDER — Final benchmark graph**
>
> This will be populated after the dedicated V1 profiling milestone.

<!--
<p align="center">
  <img src="docs/assets/cpu-gpu-benchmark.png" width="85%" alt="AstraFlow CPU versus CUDA benchmark">
</p>
-->

---

# Why I Built AstraFlow

AstraFlow began as an attempt to build a project at the intersection of several areas I wanted to understand more deeply:

- numerical mathematics;
- fluid mechanics;
- partial differential equations;
- rocket propulsion;
- modern C++;
- GPU architecture;
- CUDA programming;
- scientific software;
- verification and numerical analysis.

Rather than building a CFD visualisation around an existing solver, my goal was to implement the numerical machinery itself.

That distinction is central to the project.

AstraFlow does **not** call an external CFD package to obtain a solution and then display it. The finite-volume solver, flux calculations, reconstruction, time integration, axisymmetric geometry, transport terms, CUDA execution path, engineering integrations, verification problems, and analysis infrastructure are implemented as part of AstraFlow.

The resulting project is intended not only to produce fluid-flow images, but also to answer a deeper engineering question:

> **Can a compact, independently implemented CFD system accurately reproduce the fundamental physics of compressible rocket-nozzle flow while making effective use of a modern consumer NVIDIA GPU?**

---

# Development Philosophy

AstraFlow is being developed using a deliberately verification-first process.

The development order is approximately:

```text
Mathematical model
        ↓
1D reference problem
        ↓
CPU finite-volume solver
        ↓
Analytical verification
        ↓
CUDA implementation
        ↓
CPU ↔ GPU numerical parity
        ↓
2D conservative finite volumes
        ↓
Axisymmetric formulation
        ↓
Rocket-nozzle geometry
        ↓
Viscous Navier–Stokes transport
        ↓
Engineering quantities
        ↓
Interactive visualisation
        ↓
Profiling and optimisation
        ↓
Final scientific verification
```

The project therefore does not accept:

```text
"the simulation looks correct"
```

as evidence that the implementation is correct.

Each major numerical subsystem is tested against either:

- a known analytical solution;
- a canonical CFD verification problem;
- conservation properties;
- the independently executable CPU backend;
- mesh-refinement behaviour.

Only after these gates pass is functionality promoted into the main application.

---

# Human + AI Development Workflow

AstraFlow is also an experiment in using modern AI systems effectively during the development of technically demanding scientific software.

The project is **AI-assisted rather than AI-autonomous**.

I act as the project's architect and mathematical/engineering decision-maker. My responsibilities include:

- defining the physical scope;
- selecting the governing equations;
- selecting the numerical formulation;
- deriving and reviewing mathematical relationships;
- determining what should and should not be included in V1;
- defining verification problems;
- evaluating numerical results;
- specifying CPU/GPU parity requirements;
- deciding the CUDA architecture;
- defining milestone acceptance criteria;
- reviewing implementation decisions;
- identifying scientifically misleading shortcuts;
- directing iterative refinement.

AI coding agents are used as implementation collaborators for tasks such as:

- translating mathematical formulations into C++/CUDA;
- generating initial implementations;
- developing tests;
- investigating compiler errors;
- implementing repetitive infrastructure;
- running verification cases;
- analysing failures;
- performing targeted refactoring;
- documenting implemented systems.

The workflow is deliberately not:

```text
prompt → generate repository → assume it works
```

Instead:

```text
Architecture
    ↓
Mathematics
    ↓
Implementation specification
    ↓
AI-assisted implementation
    ↓
Numerical tests
    ↓
Human/architectural review
    ↓
Failure analysis
    ↓
Targeted correction
    ↓
Verification
    ↓
Next milestone
```

For a project involving nonlinear PDEs and GPU numerical computing, this distinction matters.

An implementation can compile successfully and produce visually plausible flow fields while still containing serious errors in:

- conservation;
- boundary conditions;
- geometry terms;
- Riemann fluxes;
- floating-point behaviour;
- dimensional scaling;
- CUDA execution;
- numerical stability.

AstraFlow therefore treats AI-generated code as something that must be **tested against mathematics and physical invariants**, not as an automatically authoritative result.

---

# Scientific Scope

The V1 solver models a:

> **single-species, calorically perfect, compressible ideal gas flowing through an axisymmetric converging-diverging rocket nozzle.**

The current solver includes:

- compressible flow;
- conservation of mass;
- axial momentum;
- radial momentum;
- total energy;
- pressure-density-energy coupling;
- finite-volume spatial discretisation;
- second-order MUSCL reconstruction;
- HLLC approximate Riemann fluxes;
- HLLE robustness fallback;
- explicit SSP-RK2 time integration;
- CFL-controlled timesteps;
- axisymmetric geometric terms;
- Newtonian viscosity;
- Fourier heat conduction;
- configurable thermodynamic properties;
- slip and no-slip wall treatment;
- stagnation inlet conditions;
- pressure/supersonic outlet handling.

The V1 model intentionally excludes:

- combustion chemistry;
- multi-species reacting flow;
- turbulence modelling;
- LES/DNS turbulence research;
- external rocket plumes;
- adaptive mesh refinement;
- conjugate wall heat transfer;
- regenerative cooling;
- ablative materials;
- non-equilibrium chemistry;
- real-gas thermodynamics;
- multi-GPU execution.

These are possible future extensions rather than shortcuts to be inserted before the baseline solver has been fully validated.

---

# Governing Equations

AstraFlow numerically solves the conservation equations for a compressible fluid.

The conservative state is represented by

\[
\mathbf U =
\begin{bmatrix}
\rho \\
\rho u \\
\rho v \\
\rho E
\end{bmatrix},
\]

where:

- \(\rho\) is density;
- \(u\) is axial velocity;
- \(v\) is radial velocity;
- \(E\) is total specific energy.

At a high level, the conservation law can be written as

\[
\frac{\partial \mathbf U}{\partial t}
+
\nabla\cdot \mathbf F_c
-
\nabla\cdot \mathbf F_v
=
\mathbf S,
\]

where:

- \(\mathbf F_c\) represents convective fluxes;
- \(\mathbf F_v\) represents viscous and thermal transport;
- \(\mathbf S\) contains the terms required by the axisymmetric formulation.

For the inviscid Euler limit,

\[
\mathbf F_v=0.
\]

For the Navier–Stokes formulation, viscous stresses and conductive heat flux are included.

---

## Continuity

Mass conservation is governed by

\[
\frac{\partial \rho}{\partial t}
+
\nabla\cdot(\rho\mathbf u)
=
0.
\]

---

## Momentum

The compressible momentum equation is

\[
\frac{\partial(\rho\mathbf u)}{\partial t}
+
\nabla\cdot
\left(
\rho\mathbf u\otimes\mathbf u+p\mathbf I-\boldsymbol\tau
\right)
=
\mathbf S_m.
\]

Here:

- \(p\) is static pressure;
- \(\boldsymbol\tau\) is the viscous stress tensor;
- \(\mathbf S_m\) contains the appropriate axisymmetric geometric contribution.

---

## Energy

Total energy evolves according to

\[
\frac{\partial(\rho E)}{\partial t}
+
\nabla\cdot
\left[
(\rho E+p)\mathbf u
-
\boldsymbol\tau\cdot\mathbf u
+
\mathbf q
\right]
=
S_E.
\]

Conductive heat transfer follows Fourier's law,

\[
\mathbf q=-k\nabla T.
\]

---

## Equation of State

V1 assumes a calorically perfect ideal gas,

\[
p=(\gamma-1)\rho e,
\]

where \(e\) is internal energy.

The local speed of sound is

\[
a=\sqrt{\frac{\gamma p}{\rho}}
\]

and Mach number is

\[
M=\frac{|\mathbf u|}{a}.
\]

This approximation is intentionally simpler than the thermochemical behaviour of real rocket exhaust.

The objective of V1 is to establish a trustworthy CFD and CUDA foundation before introducing substantially more complex gas models.

---

# Why Navier–Stokes?

Navier–Stokes equations describe conservation of momentum in viscous fluids and form the basis of an enormous amount of computational fluid dynamics.

AstraFlow **numerically approximates a particular compressible Navier–Stokes problem**.

It does not claim to solve the mathematical Navier–Stokes existence-and-smoothness problem.

The famous three-dimensional existence/smoothness question is a fundamentally different mathematical problem from numerically solving discretised Navier–Stokes equations for specified geometries, initial conditions, constitutive assumptions, and boundary conditions.

This distinction is maintained explicitly throughout the project.

---

# Numerical Method

## Finite-Volume Formulation

AstraFlow uses a **cell-centred conservative finite-volume method**.

For a control volume \(V_i\),

\[
\frac{d\mathbf U_i}{dt}
=
-
\frac{1}{V_i}
\sum_f
\mathbf F_fA_f
+
\mathbf S_i.
\]

The finite-volume formulation was selected because conservation laws are fundamental to compressible fluid dynamics.

It is particularly suitable for problems containing:

- shocks;
- rarefaction waves;
- contact discontinuities;
- strong pressure gradients;
- transonic acceleration.

Rocket nozzles can contain several of these behaviours simultaneously.

---

# MUSCL Reconstruction

First-order finite-volume methods are robust but excessively diffusive.

AstraFlow therefore implements **MUSCL reconstruction** to obtain higher spatial accuracy while retaining control around discontinuities.

The solver provides slope limiting strategies including:

- Minmod;
- Van Leer;
- Monotonized Central / MC.

These limiters reduce nonphysical oscillations around shocks and strong gradients.

---

# HLLC and HLLE Fluxes

Fluxes between neighbouring finite-volume cells require solving, exactly or approximately, a local Riemann problem.

AstraFlow uses the **HLLC approximate Riemann solver** as the primary convective flux.

HLLC is especially useful for compressible-flow problems because it resolves the approximate wave structure associated with:

- left-moving waves;
- contact waves;
- right-moving waves.

A more diffusive **HLLE fallback** is available when a reconstructed HLLC state becomes numerically pathological.

The fallback is treated as a numerical safety mechanism rather than the default solver.

Usage can be monitored during simulation so that a case requiring excessive fallback behaviour can be identified rather than silently accepted.

---

# Time Integration

AstraFlow uses explicit **second-order Strong Stability Preserving Runge–Kutta integration (SSP-RK2)**.

The timestep is limited according to a CFL condition related to local signal propagation speeds.

Conceptually,

\[
\Delta t
\propto
\mathrm{CFL}
\frac{\Delta x}
{|\mathbf u|+a}.
\]

Viscous simulations additionally account for the more restrictive diffusive stability requirement where necessary.

---

# Axisymmetric Rocket Geometry

A full three-dimensional nozzle would be unnecessarily expensive for an initial rotationally symmetric rocket-nozzle study.

AstraFlow therefore models the meridional \(x-r\) plane.

Rotating this domain around the centreline represents the corresponding axisymmetric 3D nozzle.

Conceptually:

```text
radius r
   ↑
   │
   │ chamber
   │───────────────────╲
   │                    ╲
   │                     ╲
   │                      ╲ throat
   │                       ╲
   │                        ╲
   │                         ╲────────────── exit
   │
   └────────────────────────────────────────→ x
```

The nozzle is parameterised using quantities such as:

- chamber radius;
- chamber length;
- throat radius;
- contraction length;
- exit radius;
- expansion length;
- axial grid resolution;
- radial grid resolution.

The mesh follows the nozzle geometry rather than pretending that the flow occurs inside a rectangular Cartesian box.

---

# Boundary Conditions

V1 supports the boundary conditions required for internal nozzle simulation.

### Centreline

Axisymmetric symmetry is enforced along

\[
r=0.
\]

Special treatment prevents singular behaviour from naïvely evaluating geometric terms containing \(1/r\).

### Nozzle Wall

The inviscid solver can use a slip wall.

The viscous solver supports an:

- impermeable;
- no-slip;
- adiabatic

wall.

### Chamber / Inlet

The inlet is defined through reservoir/stagnation quantities such as:

\[
P_0,\qquad T_0.
\]

### Outlet

The outlet treatment distinguishes between subsonic and supersonic behaviour.

A prescribed ambient/back pressure can influence subsonic information entering the domain, while a correctly supersonic exit should not have every primitive quantity artificially imposed.

---

# CPU Reference Backend

AstraFlow intentionally contains a CPU solver even though GPU execution is one of the main project objectives.

The CPU backend acts as a **numerical reference implementation**.

Its purpose is not maximum CPU performance.

Its purpose is to provide:

- independent execution;
- easier debugging;
- deterministic test cases;
- verification;
- CPU/CUDA numerical comparison;
- CPU-only continuous integration.

This allows CUDA kernels to be checked against the same physical problem instead of being trusted merely because they execute successfully.

---

# CUDA Backend

The production solver is implemented using **CUDA C++**.

The primary development GPU is:

```text
NVIDIA GeForce RTX 5070 Laptop GPU
Architecture: NVIDIA Blackwell
Compute capability: 12.0
Target: sm_120
```

The development environment uses a project-local CUDA 12.8 toolchain capable of generating native Blackwell `sm_120` code.

AstraFlow's initial CUDA device test does more than query a driver property: it compiles and executes a CUDA kernel and reports the architecture used by the running binary.

This protects against a configuration where CUDA appears available but the solver is not actually being built for the intended GPU generation.

---

# GPU Memory Architecture

GPU field data is designed around **Structure of Arrays (SoA)** storage.

Conceptually:

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

This layout is chosen to favour coalesced GPU memory access when neighbouring threads access the same physical quantity for adjacent cells.

Device resources are wrapped with explicit ownership rather than scattering uncontrolled `cudaMalloc` and `cudaFree` calls throughout the codebase.

Reusable buffers are allocated outside the timestep loop wherever practical.

---

# Precision

AstraFlow supports both single- and double-precision numerical paths where appropriate.

### FP32

Primary interactive/performance mode.

This is particularly suitable for a consumer RTX GPU where FP32 throughput is substantially more important than high-rate FP64 performance.

### FP64

Verification/reference GPU mode.

It is used to investigate numerical parity between the CPU and CUDA implementations.

The solver also applies dimensional scaling/non-dimensionalisation to improve numerical conditioning rather than assuming raw SI quantities will always behave well under FP32 arithmetic.

---

# Current Verification Results

The following results come from the actual development milestone tests completed so far.

They are **measured development results**, not theoretical targets.

---

## Sod Shock Tube

The canonical Sod shock tube exercises:

- a shock wave;
- a contact discontinuity;
- a rarefaction fan.

For a 400-cell CPU case, the current implementation produced mean absolute errors of approximately:

| Quantity | Mean absolute error |
|---|---:|
| Density | `0.001446` |
| Velocity | `0.002546` |
| Pressure | `0.000891` |

Mass conservation error was approximately:

```text
3.33 × 10^-15
```

for the tested configuration.

---

## CPU ↔ CUDA Euler Parity

After matching floating-point contraction behaviour between CPU and CUDA code paths:

### FP64

Maximum CPU/GPU difference:

```text
≈ 4.22 × 10^-14
```

### FP32

Maximum CPU/GPU difference:

```text
≈ 4.62 × 10^-6
```

Both precision modes passed the configured verification tolerances.

---

## 2D Conservation Test

A transported 2D density-wave case conserved all four integrated conservative quantities within approximately:

```text
1 × 10^-12
```

Maximum CPU/GPU difference was approximately:

```text
1.78 × 10^-15
```

for the tested FP64 configuration.

---

## Axisymmetric Nozzle Verification

A smooth inviscid nozzle verification case currently produces:

```text
Throat Mach number: ≈ 1.0096
```

showing the expected transition through approximately sonic conditions at the throat.

The mean Mach-number difference from the corresponding quasi-one-dimensional isentropic relation was approximately:

```text
0.000319
```

and axial mass-flow spread was approximately:

```text
0.0774 %
```

for the tested configuration.

A stationary axisymmetric state was preserved with a maximum measured drift of approximately:

```text
2.22 × 10^-15
```

and nozzle CPU/GPU parity differed by approximately:

```text
4.88 × 10^-15
```

in the tested configuration.

---

## Viscous Verification

Compressible Couette flow is used as a canonical viscous verification problem.

Current mean errors are approximately:

| Quantity | Mean error |
|---|---:|
| Velocity | `5.15 × 10^-5` |
| Temperature | `4.46 × 10^-5` |

CPU/GPU viscous-state parity reached approximately:

```text
FP64: 4.00 × 10^-15
FP32: 4.26 × 10^-6
```

for the milestone validation cases.

---

# Verification Still Required Before V1 Acceptance

These development results are encouraging, but V1 is not considered scientifically complete yet.

The final verification phase will additionally include:

- complete Release CPU test suite;
- complete Release CUDA test suite;
- final uniform-state verification;
- Sod regression;
- isentropic nozzle regression;
- viscous analytical regression;
- CPU/GPU parity regression;
- multi-resolution grid-refinement study;
- rocket-nozzle runtime tests;
- CUDA memory/error diagnostics;
- final engineering-output consistency analysis.

AstraFlow will not claim a formal spatial convergence result until the planned grid-refinement study has actually been completed and documented.

---

# Rocket Engineering Outputs

AstraFlow transforms numerical flow fields into useful engineering measurements.

## Mass Flow

Mass flow is evaluated by integrating the numerical field across an axisymmetric section,

\[
\dot m
=
\int_A \rho u\,dA.
\]

Rather than assuming the exit is perfectly uniform, the integration uses the finite-volume solution.

---

## Exit Conditions

The solver determines quantities including:

- exit Mach number;
- axial velocity;
- pressure;
- temperature.

Appropriate area- or mass-weighted averaging is used depending on the quantity.

---

## Thrust

Nozzle thrust is estimated directly from the computed exit plane,

\[
F
=
\int_{A_e}\rho u^2\,dA
+
\int_{A_e}(p-p_a)\,dA.
\]

This allows nonuniform exit profiles to contribute naturally.

---

## Specific Impulse

Specific impulse is calculated as

\[
I_{sp}
=
\frac{F}{\dot m g_0}.
\]

These values describe the idealised CFD model and must not be interpreted as experimental certification of a real rocket engine.

---

# Reproducible Simulation Output

AstraFlow is designed so that simulations can be reproduced and inspected independently of the GUI.

A run can produce:

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

Stores the effective simulation configuration.

### `convergence.csv`

Stores residual/convergence history.

### `performance.json`

Stores measured runtime information.

### `summary.json`

Stores engineering outputs and simulation metadata.

### `mesh.vts`

Exports the structured nozzle mesh.

### `final_state.vts`

Exports CFD fields for independent visualisation.

VTK StructuredGrid output allows the numerical results to be loaded into tools such as **ParaView**, providing an independent visualisation path instead of forcing all analysis through AstraFlow's own interface.

---

# Interactive GUI

The V1 desktop application uses:

- Dear ImGui;
- ImPlot;
- GLFW;
- OpenGL.

The GUI is deliberately separated from the CFD implementation.

There is no duplicate GUI-specific physics solver.

The interface invokes the same underlying solver infrastructure used by command-line simulations.

The application is designed around panels for:

### Simulation Control

- Run;
- Pause;
- Single Step;
- Reset;
- Regenerate Mesh.

### Numerical Parameters

- CPU/CUDA backend;
- precision;
- grid resolution;
- CFL;
- iteration limit;
- residual target.

### Gas / Chamber

- \(\gamma\);
- gas constant;
- viscosity;
- Prandtl number;
- stagnation pressure;
- stagnation temperature;
- ambient pressure.

### Geometry

- chamber size;
- throat radius;
- contraction dimensions;
- exit radius;
- expansion length.

### Flow Visualisation

Fields include:

- pressure;
- density;
- temperature;
- Mach number;
- axial velocity;
- radial velocity;
- velocity magnitude;
- total energy;
- vorticity.

### Convergence

Live residual plots display the evolution of:

- continuity;
- axial momentum;
- radial momentum;
- energy.

### Engineering Analysis

The interface exposes:

- mass flow;
- throat Mach;
- exit Mach;
- exit velocity;
- exit pressure;
- thrust;
- specific impulse;
- conservation diagnostics.

### GPU Performance

The final interface is designed to display:

- detected NVIDIA GPU;
- CUDA architecture;
- grid cell count;
- device memory;
- iteration time;
- iterations per second;
- visualisation update rate.

---

# Asynchronous Simulation Architecture

The GUI is not intended to block while CFD iterations execute.

Simulation execution and user-interface rendering are therefore separated.

The solver exposes states conceptually equivalent to:

```text
Idle
Ready
Running
Paused
Converged
Failed
```

The GUI communicates with the solver worker and receives safe snapshots of simulation fields.

This prevents the rendering thread from arbitrarily modifying CUDA state while kernels are operating.

For V1, selected visualisation data can be periodically transferred from device memory to host memory.

Direct CUDA/OpenGL interoperability is intentionally deferred until profiling demonstrates that these transfers are a meaningful bottleneck.

---

# Command-Line Interface

AstraFlow also provides a headless CLI.

Conceptually:

```bash
astraflow_cli \
    --config examples/rocket_nozzle/config.json \
    --backend cuda
```

Backends include:

```text
cpu
cuda
```

A command-line simulation can report:

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
Estimated thrust
Specific impulse
Mass conservation error
Final residual
```

This separation is important because scientific verification and performance benchmarking should not require a graphical desktop.

---

# Performance Objectives

A major objective of AstraFlow is to investigate where GPU execution becomes advantageous for finite-volume CFD.

The benchmark programme will compare representative grid sizes such as:

```text
128 × 32
256 × 64
512 × 128
1024 × 256
```

where practical.

Measurements will include:

| Metric | Purpose |
|---|---|
| Grid cells | problem scale |
| CPU iteration time | reference performance |
| CUDA iteration time | GPU performance |
| Iterations/s | solver throughput |
| CPU/GPU speedup | acceleration |
| GPU memory | memory scaling |
| Kernel time | optimisation target |

No arbitrary minimum speedup is being imposed.

For small meshes, CUDA may legitimately be slower than the CPU because launch, synchronisation, and data-management overhead dominate the amount of useful arithmetic.

The expected crossover toward GPU advantage should occur as the problem becomes sufficiently large and parallel.

The final benchmark report will therefore show the **measured scaling behaviour**, even if it is less impressive than originally expected.

---

# CUDA Optimisation Strategy

Optimisation is performed only after correctness.

The planned workflow is:

```text
Measure
   ↓
Profile
   ↓
Identify dominant kernel / memory cost
   ↓
Optimise targeted bottleneck
   ↓
Rerun CPU ↔ GPU parity
   ↓
Rerun scientific verification
   ↓
Measure again
```

Likely areas of investigation include:

- memory coalescing;
- global-memory traffic;
- temporary buffers;
- kernel-launch overhead;
- reduction kernels;
- warp divergence;
- occupancy;
- unnecessary synchronisation;
- host/device transfers.

Where available, profiling will use NVIDIA tooling such as Nsight Systems, Nsight Compute, and Compute Sanitizer.

Performance changes are not accepted if they invalidate numerical verification.

---

# Expected V1 Results

The final V1 should be capable of demonstrating several characteristic features of compressible nozzle flow.

These are **expected physical behaviours**, not pre-recorded result claims.

## Choking

For sufficiently high chamber-to-ambient pressure ratios, the flow should accelerate toward approximately

\[
M=1
\]

around the nozzle throat.

This behaviour is already visible in the current inviscid nozzle verification case.

---

## Supersonic Expansion

After the throat, the diverging section should support acceleration to

\[
M>1
\]

under suitable pressure conditions.

The final GUI should make this immediately visible through Mach-number contours.

---

## Pressure Conversion

The nozzle should demonstrate conversion of:

```text
high pressure / thermal energy
        ↓
directed kinetic energy
```

as gas moves from the chamber through the converging-diverging geometry.

---

## Back-Pressure Sensitivity

Changing ambient/back pressure should alter the nozzle solution.

The project is intended eventually to illustrate conditions corresponding to:

- approximately ideal expansion;
- over-expansion;
- under-expansion;
- internal compression/shock behaviour where captured by the model and computational domain.

---

## Viscous Effects

Compared with the inviscid solution, the Navier–Stokes model should reveal effects associated with:

- near-wall velocity gradients;
- viscous stress;
- thermal conduction;
- boundary-layer behaviour;
- changes in effective nozzle performance.

---

# Project Architecture

The repository follows a modular layout conceptually similar to:

```text
AstraFlow/
│
├── apps/
│   ├── cli/
│   └── gui/
│
├── include/
│   └── astraflow/
│
├── src/
│   ├── core/
│   ├── cpu/
│   ├── cuda/
│   │   ├── kernels/
│   │   ├── memory/
│   │   └── reductions/
│   ├── physics/
│   ├── numerics/
│   ├── geometry/
│   ├── analysis/
│   ├── io/
│   └── visualization/
│
├── tests/
│   ├── unit/
│   ├── verification/
│   └── regression/
│
├── benchmarks/
│
├── examples/
│
├── scripts/
│
├── docs/
│   ├── mathematics/
│   ├── architecture/
│   ├── verification/
│   └── performance/
│
└── .github/
    └── workflows/
```

The design maintains separation between:

```text
Physics
Numerics
Geometry
CPU backend
CUDA backend
Analysis
I/O
Visualisation
Applications
```

One of the project's architectural rules is:

> **CFD mathematics does not belong inside GUI code.**

---

# Technology Stack

## Core

- C++20
- CMake
- Ninja

## GPU

- NVIDIA CUDA C++
- native Blackwell `sm_120` target
- CUDA 12.8-class toolchain

## Interface

- Dear ImGui
- ImPlot
- GLFW
- OpenGL

## Testing

- Catch2

## Configuration

- nlohmann/json

## Scientific Output

- JSON
- CSV
- VTK XML StructuredGrid (`.vts`)

## Auxiliary Tooling

Python may be used for:

- reference calculations;
- plotting;
- benchmark processing;
- scientific analysis.

It is not the production CFD implementation.

---

# Build

## Requirements

The CPU build requires approximately:

```text
CMake >= 3.24
Ninja
C++20-compatible compiler
```

The CUDA build additionally requires a CUDA toolchain capable of targeting the GPU architecture selected at configuration time.

For the primary AstraFlow development machine:

```text
NVIDIA RTX 5070 Laptop GPU
Blackwell
Compute Capability 12.0
sm_120
```

---

## CPU Build

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

## CUDA Build

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

## Run Tests

```bash
ctest \
    --test-dir build \
    --output-on-failure
```

CPU-only verification can similarly use:

```bash
ctest \
    --test-dir build-cpu \
    --output-on-failure
```

---

# Dependencies

The dependency policy is intentionally conservative.

Current lightweight dependencies include:

| Dependency | Purpose |
|---|---|
| `nlohmann/json 3.11.3` | simulation configuration and structured output |
| `Catch2 3.7.1` | unit and verification testing |
| Dear ImGui | interactive interface |
| ImPlot | numerical/convergence plots |
| GLFW | window/input handling |
| OpenGL | visualisation rendering |
| CUDA Runtime | GPU execution |

Upstream licenses remain applicable to their respective dependencies.

AstraFlow itself is released under the MIT License.

---

# Development Milestones

## Milestone 0 — Toolchain and Repository

**Status: COMPLETE**

- project scaffold;
- CMake;
- CPU compilation;
- CUDA compilation;
- RTX 5070 execution;
- native `sm_120` validation;
- test infrastructure.

---

## Milestone 1 — CPU Euler Solver

**Status: COMPLETE**

- primitive/conservative state;
- thermodynamics;
- HLLC;
- HLLE fallback;
- MUSCL;
- limiters;
- CFL;
- SSP-RK2;
- Sod shock tube;
- conservation tests.

---

## Milestone 2 — CUDA Euler Solver

**Status: COMPLETE**

- GPU timestep execution;
- FP32;
- FP64;
- CPU/GPU parity;
- CUDA validation.

---

## Milestone 3 — 2D Finite Volumes

**Status: COMPLETE**

- 2D cells/faces;
- gradients;
- multidimensional fluxes;
- CPU backend;
- CUDA backend;
- conservation verification.

---

## Milestone 4 — Axisymmetric Rocket Nozzle

**Status: COMPLETE**

- nozzle geometry;
- body-conforming structured mesh;
- axisymmetric terms;
- inlet/outlet treatment;
- sonic throat behaviour;
- quasi-1D comparison.

---

## Milestone 5 — Compressible Navier–Stokes

**Status: COMPLETE**

- Newtonian viscous stress;
- gradients;
- viscosity;
- heat conduction;
- no-slip adiabatic wall;
- Couette verification;
- CUDA parity.

---

## Milestone 6 — Engineering Analysis / CLI / Output

**Status: COMPLETE**

- configuration loading;
- CLI;
- mass flow;
- exit quantities;
- thrust;
- specific impulse;
- conservation metrics;
- residual history;
- JSON;
- CSV;
- VTK StructuredGrid.

---

## Milestone 7 — Interactive GUI

**Status: IN PROGRESS**

The implementation includes the solver worker/state architecture and GUI controls. Final rendering and WSLg acceptance testing remain underway.

---

## Milestone 8 — CUDA Performance Pass

**Status: PENDING**

Will include:

- Release benchmarking;
- CUDA profiling;
- bottleneck identification;
- targeted optimisation;
- post-optimisation scientific regression.

---

## Milestone 9 — Final Phase 1 Acceptance

**Status: PENDING**

Will include:

- complete CPU build;
- complete CUDA build;
- verification suite;
- grid refinement;
- CPU/GPU rocket runs;
- GUI acceptance;
- runtime diagnostics;
- documentation;
- repository audit.

---

# Continuous Integration

AstraFlow includes CPU-focused CI so that the numerical core can be compiled and tested without requiring an NVIDIA GPU on the CI runner.

The CI pipeline is intended to cover:

```text
CMake configuration
CPU compilation
unit tests
small verification tests
```

CUDA verification and GPU performance testing are performed on the actual Blackwell development hardware.

---

# Scientific Integrity

AstraFlow is a portfolio project, but numerical credibility takes precedence over producing impressive screenshots or benchmark claims.

The project therefore follows several rules.

### No fabricated validation

A test is not marked as passing merely because a result looks plausible.

### No fabricated acceleration

CUDA speedup figures are reported from actual measurements.

### No artificially loose tolerances

Numerical tolerances should represent reasonable floating-point and discretisation behaviour rather than values chosen simply to make CI green.

### No hidden instability

NaN, negative-density, negative-pressure, and invalid-state behaviour is explicitly detected.

### No overstated physical fidelity

A calorically perfect ideal-gas nozzle simulation is not presented as a complete simulation of a chemically reacting rocket engine.

### No claim of solving the Navier–Stokes Millennium Problem

Numerical CFD and the mathematical existence/smoothness question are distinct problems.

---

# Limitations of V1

The first AstraFlow release should be understood as a **verified CFD foundation**, not as a production aerospace design package.

Important limitations include:

- axisymmetric rather than general 3D geometry;
- ideal-gas thermodynamics;
- single species;
- no combustion chemistry;
- no turbulence closure;
- no external exhaust plume;
- no adaptive mesh refinement;
- no wall-material heat conduction;
- no experimental calibration;
- no certification-grade validation;
- consumer-GPU-oriented architecture.

These constraints are deliberate.

They create a sufficiently complex system for serious CFD/HPC work while retaining a verification scope that can be independently understood and tested.

---

# Future Development

Once V1 has passed its final acceptance review, possible future directions include:

## External Plume Simulation

Extend the domain beyond the nozzle exit to investigate:

- under-expanded jets;
- over-expanded jets;
- expansion structures;
- shock cells;
- plume interaction.

## Higher-Fidelity Thermodynamics

Introduce:

- temperature-dependent specific heats;
- improved viscosity models;
- species thermodynamics.

## Reacting Flow

Potentially incorporate:

- multiple species;
- reaction kinetics;
- combustion modelling.

## Turbulence

Investigate models such as:

- RANS;
- \(k-\omega\) SST;
- eventually LES for suitable research cases.

## Thermal Analysis

Extend toward:

- wall heat flux;
- conjugate heat transfer;
- regenerative cooling studies.

## Advanced GPU Optimisation

Potential directions include:

- CUDA/OpenGL interoperability;
- more aggressive kernel fusion;
- asynchronous execution;
- improved reduction algorithms;
- memory-layout experiments;
- multi-GPU execution.

## Nozzle Optimisation

Use the CFD solver as the evaluation engine for parameterised nozzle-design studies.

Potential objectives include:

- thrust;
- specific impulse;
- pressure matching;
- wall loading;
- efficiency.

---

# What I Am Learning Through AstraFlow

AstraFlow is intentionally broader than a single programming exercise.

It provides practical experience with:

### Applied Mathematics

- conservation laws;
- PDE discretisation;
- Riemann problems;
- stability;
- truncation error;
- convergence;
- dimensional analysis.

### Fluid Dynamics

- compressible flow;
- shocks;
- sonic conditions;
- nozzle expansion;
- viscosity;
- heat conduction.

### Numerical Methods

- finite-volume schemes;
- reconstruction;
- slope limiting;
- Runge–Kutta integration;
- CFL stability;
- verification methodology.

### High-Performance Computing

- CUDA kernels;
- GPU memory;
- floating-point behaviour;
- parallel reductions;
- profiling;
- CPU/GPU scaling.

### Software Engineering

- modular C++;
- CMake;
- backend abstraction;
- configuration;
- automated testing;
- reproducible output;
- CI;
- technical documentation.

### AI-Assisted Engineering

Perhaps most importantly, the project is helping me investigate how AI coding systems can be used productively in technical work without outsourcing understanding.

The objective is not to remove myself from the mathematical or engineering process.

It is to use AI to accelerate implementation while retaining responsibility for:

- the model;
- the assumptions;
- the architecture;
- the tests;
- the interpretation;
- the final technical judgement.

---

# Repository Status

Phase 1 development is ongoing.

For the detailed live milestone ledger, see:

[`docs/PHASE1_STATUS.md`](docs/PHASE1_STATUS.md)

The README will be updated with:

- final GUI screenshots;
- final CFD visualisations;
- complete grid-convergence data;
- CPU/GPU benchmark plots;
- profiling results;
- final V1 acceptance metrics

after the corresponding milestones have passed.

---

# License

AstraFlow is licensed under the **MIT License**.

Third-party software retains its respective upstream licensing.

NVIDIA CUDA and associated NVIDIA tooling remain subject to NVIDIA's applicable license terms.

---

# Disclaimer

AstraFlow is an educational, research, and portfolio scientific-computing project.

It is **not** flight-qualified aerospace software and should not be used as the sole basis for safety-critical rocket-engine design, manufacturing, or operation.

Its results should be interpreted according to the assumptions, discretisation, numerical error, verification evidence, and physical limitations documented in the repository.

---

## Author

**Arnav**

BSc Computer Science with Artificial Intelligence  
University of Nottingham

Interests:

- scientific computing;
- artificial intelligence;
- numerical methods;
- GPU computing;
- high-performance computing;
- computational science.

---

> **AstraFlow's objective is not simply to make fluid move on a screen.**
>
> It is to understand how the mathematics of conservation laws becomes reliable numerical software, how that software maps onto massively parallel GPU hardware, and how the resulting computation can be verified before its output is trusted.

# Architecture evolution

This document describes the intended architectural direction for future AstraFlow releases. It supplements [overview.md](overview.md), which remains the source for the **current implemented architecture**.

## Current boundary

The V1 design has a useful separation:

~~~mermaid
flowchart LR
    M["Physics + numerics"] --> B["CPU / CUDA backends"]
    B --> S["Simulation"]
    S --> A["Analysis / I/O"]
    S --> U["CLI / GUI"]
~~~

That separation should survive Phase 2. The largest future risk is allowing new physics to appear as backend-specific conditionals or GUI-owned calculations.

## Target subsystem boundaries

~~~mermaid
flowchart TB
    subgraph Models
        TH["Thermodynamics"]
        TR["Transport"]
        TURB["Turbulence (later)"]
        CHEM["Species/Chemistry (later)"]
    end

    subgraph Discretization
        GEO["Geometry / mesh"]
        REC["Reconstruction"]
        FLX["Convective flux"]
        VIS["Viscous flux"]
        SRC["Axisymmetric/source terms"]
        BC["Boundary models"]
        RES["Residual assembly R(U)"]
    end

    subgraph Execution
        CPU["CPU scheduler"]
        CUDA["CUDA scheduler"]
        STEADY["Steady accelerator"]
        TRANS["Physical-time integrator"]
    end

    subgraph Products
        ANA["Engineering / wall / plume analysis"]
        IO["Run artifacts"]
        GUI["Immutable GUI snapshots"]
    end

    TH --> FLX
    TH --> VIS
    TR --> VIS
    TURB -. future .-> VIS
    CHEM -. future .-> RES

    GEO --> RES
    REC --> FLX
    FLX --> RES
    VIS --> RES
    SRC --> RES
    BC --> FLX

    RES --> CPU
    RES --> CUDA
    CPU --> TRANS
    CUDA --> TRANS
    CPU --> STEADY
    CUDA --> STEADY

    TRANS --> ANA
    STEADY --> ANA
    ANA --> IO
    ANA --> GUI
~~~

The architectural invariant is that **physical-time and steady acceleration evaluate the same spatial residual for the same model/configuration**. Their difference is how they march toward a state, not the equations they claim to solve.

## Model interfaces

Future thermodynamics/transport abstraction should satisfy three competing constraints:

1. readable scientific API;
2. CPU/CUDA parity;
3. no expensive dynamic dispatch in hot device code.

Prefer plain data/configuration plus compile-time or enum-selected device-callable functions over deep polymorphic hierarchies.

Conceptually:

~~~text
GasModelData
  kind
  constants / table coefficients

thermo_from_conservative(U, model)
thermo_from_T_rho(T, rho, model)

TransportModelData
  kind
  coefficients

mu(T, model)
k(T, thermo, model)
~~~

The exact implementation should follow profiling and CUDA compilation constraints.

## Mesh evolution

### V1

Single structured axisymmetric body-conforming nozzle domain.

### Phase 2 target

A structured or multi-block structured domain that can represent both the internal nozzle and external ambient region.

~~~mermaid
flowchart LR
    C["Chamber block"] --> N["Nozzle block"]
    N --> P["Near-plume block"]
    P --> F["Far-field block"]

    C -. conservative interface .- N
    N -. conservative interface .- P
    P -. conservative interface .- F
~~~

If multi-block interfaces are added, they must have:

- one authoritative face geometry;
- equal/opposite conservative exchange;
- explicit ownership;
- parity tests;
- uniform-state preservation.

Do not introduce a fully unstructured mesh until a concrete Phase 2 requirement cannot be represented robustly by the staged design.

## Boundary architecture

Boundary conditions should become explicit model objects/data rather than problem-name branches scattered through the operator.

Desired conceptual categories:

~~~text
Axis
SlipWall
NoSlipWall
ReservoirInlet
PressureOutlet
SupersonicOutlet
FarFieldCharacteristic
BlockInterface
~~~

A boundary implementation owns **state/flux semantics**, while the finite-volume operator owns conservation and face orientation.

## Time integration vs steady acceleration

Keep two explicit concepts:

### Physical integration

- advances physical time;
- uses globally consistent physical timestep;
- required for transients;
- existing SSP-RK2 is the reference.

### Steady acceleration

- advances pseudo-time or another nonlinear iteration;
- does not represent physical time;
- targets the same R(U)=0;
- must expose separate metadata;
- must pass root-equivalence evidence before becoming a production scientific path.

This separation prevents a common CFD failure mode: reporting accelerated pseudo-time as though it were a physical transient.

## Analysis architecture

Analysis should consume immutable physical snapshots/results rather than reach into backend memory.

Future analysis namespaces can grow by scientific responsibility:

~~~text
analysis/
  engineering        existing integral quantities
  wall               pressure / shear / heat flux
  plume              centerline / plume observables
  convergence        residual + engineering stability
  uncertainty        grid/Richardson/GCI helpers
~~~

Keep analysis definitions documented. A metric name should identify whether it is:

- a cell-center estimate;
- a reconstructed numerical face flux;
- an area/mass weighted average;
- a global conservative balance.

## Output schema evolution

Future result files should carry a schema version.

Recommended metadata:

- source commit;
- executable hash;
- effective config hash;
- backend;
- precision;
- physics-model identifiers;
- mesh identifier;
- physical vs pseudo-time mode;
- termination reason;
- convergence criteria;
- unit/scaling conventions.

Never silently change the meaning of an existing field.

## Performance boundary

CUDA optimisation should remain behind the same solver interface. Features such as CUDA/OpenGL interop, CUDA Graphs, or multi-GPU may improve execution, but they must not become dependencies of the numerical model.

## Architecture acceptance questions

Before adding a new subsystem, answer:

1. What equation/model does it add?
2. What existing limiting case must remain unchanged?
3. Where is its state owned?
4. Is it backend-neutral?
5. How is it verified independently?
6. How does it appear in run provenance?
7. Can it fail without invalidating unrelated modules?
8. What evidence permits it to be called complete?

If those questions are unclear, the subsystem is not ready to implement.

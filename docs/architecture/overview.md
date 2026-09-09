# Architecture

The reusable `astraflow` library is linked by applications and Catch2 verification. `physics/state.hpp` defines ideal-gas conversions and a four-plane SoA view. `numerics/flux.hpp` supplies limiters and normal-frame HLLC/HLLE; `euler1d_ops.hpp` supplies shared reconstruction. These small mathematical headers are instantiated for float and double on CPU and CUDA.

`src/cpu/euler1d.cpp` owns the deterministic reference loop. `src/cuda/euler1d.cu` owns an independent backend execution pipeline calling the same numerical functions. The public CUDA class uses an out-of-line PIMPL; callers do not allocate CUDA memory. Reusable RAII buffers and events live in `src/cuda/memory/resources.cuh`. There are no allocations in the GPU timestep loop.

GPU stages are named primitive/CFL conversion, reconstruction/flux, residual assembly and RK update. CUB supplies the global minimum timestep reduction with preallocated scratch. Invalid-state indices and fallback counters cross the host boundary as diagnostics. The conservative fields remain on the device throughout integration; explicit snapshots copy them to the host.

CPU-only CI disables CUDA and runs small unit/verification suites. GPU checks require the local hardware and are separate CTest entries. `geometry` owns structured meshes; `numerics/finite_volume.hpp` shares boundary reconstruction, central transport gradients, face fluxes, source and CFL evaluation. CPU and CUDA 2D backends own scheduling and memory, with a common `Solver` interface.

`Config` strictly parses known JSON keys and rejects invalid settings. `Simulation` builds and scales a case, chooses a backend and restores physical units for snapshots. CLI and GUI share this application-independent orchestration. `analysis` computes engineering integrals from physical snapshots; `io` writes VTK directly without VTK libraries.

Each run stores config.json, convergence.csv, performance.json, summary.json, mesh.vts and final_state.vts. Explicit output paths must be empty, preventing accidental overwrite. Default paths use a unique timestamp under ignored runs/. The summary distinguishes convergence, end-time and iteration-limit termination. Residuals are RMS time derivatives in nondimensional solver units; times/fields/engineering results use configuration units (SI in rocket examples).

## Interactive interface

The Dear ImGui/ImPlot/GLFW/OpenGL application has logical simulation, field, convergence and engineering/performance panels. `Controller` owns a `std::jthread` solver worker, condition variable and command state. Run/Pause/Step/Reset/Regenerate never expose active CUDA arrays to the renderer. A pause becomes visible at a completed timestep boundary. Regeneration validates configuration before replacing the simulation.

Snapshots are immutable shared objects containing mesh, physical state, gas properties, metrics and residual history. The worker publishes every ten iterations and at pause/finish/reset boundaries. The GUI derives fields from the snapshot's own gas model; older snapshots remain valid across reset. Physics and field derivatives live in the library, not GUI event handlers.

The meridional field is displayed as piecewise-constant cell colours, mirrored about the axis for nozzle cases. Radial display magnification is explicitly labelled. Mouse wheel zoom and middle-button pan change the view only. Every RK stage continues independently of drawing. No CUDA/OpenGL interoperability is used in V1.

Pinned GUI dependencies: Dear ImGui 1.91.8 (MIT), ImPlot 0.16 (MIT), GLFW 3.4 (zlib/libpng). Their licenses remain in fetched sources. Ubuntu GUI headers can be extracted locally with `bash scripts/bootstrap_gui.sh`; the existing system OpenGL/X11 runtime is reused without installation or driver changes.

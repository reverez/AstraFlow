# Architecture

The reusable `astraflow` library is linked by applications and Catch2 verification. `physics/state.hpp` defines ideal-gas conversions and a four-plane SoA view. `numerics/flux.hpp` supplies limiters and normal-frame HLLC/HLLE; `euler1d_ops.hpp` supplies shared reconstruction. These small mathematical headers are instantiated for float and double on CPU and CUDA.

`src/cpu/euler1d.cpp` owns the deterministic reference loop. `src/cuda/euler1d.cu` owns an independent backend execution pipeline calling the same numerical functions. The public CUDA class uses an out-of-line PIMPL; callers do not allocate CUDA memory. Reusable RAII buffers and events live in `src/cuda/memory/resources.cuh`. There are no allocations in the GPU timestep loop.

GPU stages are named primitive/CFL conversion, reconstruction/flux, residual assembly and RK update. CUB supplies the global minimum timestep reduction with preallocated scratch. Invalid-state indices and fallback counters cross the host boundary as diagnostics. The conservative fields remain on the device throughout integration; explicit snapshots copy them to the host.

CPU-only CI disables CUDA and runs small unit/verification suites. GPU checks require the local hardware and are separate CTest entries. Graphical and multidimensional components will be added at their specified milestone gates.

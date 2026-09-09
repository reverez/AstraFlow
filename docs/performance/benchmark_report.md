# Benchmark report

## Release 2D viscous nozzle, equal FP32

Hardware: NVIDIA GeForce RTX 5070 Laptop GPU, CC 12.0, native SM120, CUDA 12.8.93. CPU is the single-thread reference backend. Both backends use FP32, the same mesh/state, MC, CFL=0.4 and constant nondimensional viscosity 1e-4 with no-slip walls. Each grid has 20 warmup iterations followed by 50 measured iterations. Setup, final snapshots and file output are excluded. CPU timing is steady-clock time inside step(); GPU iteration timing is a CUDA-event interval spanning kernels and required host synchronization. JSON also records wall time.

| Grid | Cells | CPU ms/iter | GPU ms/iter | GPU iter/s | CPU/GPU | Device MiB |
|---|---:|---:|---:|---:|---:|---:|
| 128 x 32 | 4,096 | 1.448 | 0.539 | 1856.5 | 2.69x | 1.54 |
| 256 x 64 | 16,384 | 6.194 | 0.512 | 1954.8 | 12.11x | 6.10 |
| 512 x 128 | 65,536 | 28.222 | 0.619 | 1614.3 | 45.56x | 24.32 |
| 1024 x 256 | 262,144 | 121.539 | 1.987 | 503.2 | 61.16x | 97.11 |

These are measured sample means, not guaranteed speedups or a comparison against a tuned multicore CPU solver. Laptop clock/load variation is visible between the baseline and optimized samples. The smallest grid regressed in this sample; no universal improvement is claimed. Buffer counts exclude CUDA context/driver bookkeeping and OpenGL resources.

## Evidence-based reduction change

The instrumented 512x128 baseline identified diagnostic reductions/transfers as the largest measured interval (0.306395 ms). The old pipeline wrote four residual-square planes and reduced each separately. The final pipeline transforms the existing residual into a four-component value during one CUB reduction. It eliminates the square buffer and three reduction calls. Conservative updates and CFL selection are unchanged. A test verifies that enabling instrumentation leaves the solution identical. All scientific/parity tests passed after the change.

| Grid | Baseline GPU ms | Final GPU ms | Baseline bytes | Final bytes |
|---|---:|---:|---:|---:|
| 128 x 32 | 0.475803 | 0.538640 | 1667625 | 1619751 |
| 256 x 64 | 0.643517 | 0.511556 | 6649383 | 6400295 |
| 512 x 128 | 0.664173 | 0.619459 | 26532391 | 25496871 |
| 1024 x 256 | 2.090762 | 1.987206 | 106013223 | 101831975 |

## CUDA-event stage intervals

30 measured iterations after 20 warmup iterations on 512x128. Instrumentation adds event/host overhead, so use the uninstrumented table for iteration performance. Intervals include any queue/host gaps between their markers; these are not Nsight kernel-only measurements.

| Stage | Before ms | After ms |
|---|---:|---:|
| primitive_conversion | 0.145164 | 0.108186 |
| CFL_and_host_dt | 0.127011 | 0.111877 |
| gradients | 0.137337 | 0.113068 |
| convective_viscous_faces | 0.095516 | 0.104278 |
| residual_assembly | 0.056548 | 0.055268 |
| RK_update | 0.034870 | 0.021583 |
| diagnostic_reductions_and_transfers | 0.306395 | 0.184032 |

Nsight Systems, Nsight Compute and compute-sanitizer were absent at preflight. No Nsight or compute-sanitizer validation is claimed. Optional `numerics.profile=true` records stage intervals in the solver statistics; the benchmark averages these explicitly. Raw measurements are kept in baseline.json, optimized.json, profile-baseline.json and profile-optimized.json.

## Reproduce

```sh
cmake -S . -B build -DASTRAFLOW_BUILD_BENCHMARKS=ON
cmake --build build
./build/astraflow_benchmark --output runs/benchmark.json
./build/astraflow_benchmark --profile --nx 512 --nr 128 --iterations 30 --output runs/profile.json
```

Output files must be new. Benchmarks do not assert an expected speedup. The largest tested grid was 1024x256 (262,144 cells); maximum conservative CPU/GPU differences across the measured grids remained below 8e-6.

## Initial 1D measurement

At M2, 400-cell Sod integration measured CPU/GPU FP32 iteration times 0.0523683/0.248807 ms, and FP64 0.0629333/0.255372 ms. This small 1D case was slower on GPU. Those earlier measurements use the original 1D implementation and are not directly comparable to the 2D table.

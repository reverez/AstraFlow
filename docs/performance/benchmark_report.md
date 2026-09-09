# Benchmark report

## Initial 1D baseline (M2)

Release, RTX 5070 Laptop GPU, CUDA 12.8.93, SM120, nx=400 Sod to t=0.2. CPU wall time and CUDA-event elapsed time span integration, excluding setup and final snapshots. These are individual measurements, not statistical claims.

| Precision | CPU ms/iteration | GPU ms/iteration | CPU/GPU ratio |
|---|---:|---:|---:|
| FP32 | 0.0523683 | 0.248807 | 0.210 |
| FP64 | 0.0629333 | 0.255372 | 0.246 |

This small grid is slower on GPU. Multiple launches and host/device timestep/validity synchronization dominate such a small workload; this is a pipeline-based inference, not an Nsight measurement. No production 2D speedup is claimed. Representative grid benchmarking is reserved for M8.

#include "astraflow/core/device.hpp"
#include <cuda_runtime.h>
#include <sstream>
#include <stdexcept>
namespace astraflow {
namespace {
void check(cudaError_t e) {
    if (e != cudaSuccess)
        throw std::runtime_error(cudaGetErrorString(e));
}
__global__ void architecture_probe(int *architecture) {
#ifdef __CUDA_ARCH__
    *architecture = __CUDA_ARCH__;
#endif
}
} // namespace
std::string device_info() {
    cudaDeviceProp p{};
    check(cudaGetDeviceProperties(&p, 0));
    int runtime = 0;
    check(cudaRuntimeGetVersion(&runtime));
    int *device = nullptr;
    check(cudaMalloc(&device, sizeof(int)));
    architecture_probe<<<1, 1>>>(device);
    auto launch = cudaGetLastError();
    int architecture = 0;
    auto copy = cudaMemcpy(&architecture, device, sizeof(int), cudaMemcpyDeviceToHost);
    auto release = cudaFree(device);
    check(launch);
    check(copy);
    check(release);
    std::ostringstream out;
    out << p.name << "\nCompute capability: " << p.major << '.' << p.minor
        << "\nGlobal memory (bytes): " << p.totalGlobalMem << "\nCUDA runtime: " << runtime
        << "\nExecuted kernel architecture: sm_" << architecture / 10;
    return out.str();
}
} // namespace astraflow

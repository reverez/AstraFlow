#include "astraflow/core/device.hpp"
#include "memory/resources.cuh"
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
    cuda::Buffer<int> device(1);
    architecture_probe<<<1, 1>>>(device.data());
    AF_CUDA(cudaGetLastError());
    int architecture = device.download()[0];
    std::ostringstream out;
    out << p.name << "\nCompute capability: " << p.major << '.' << p.minor
        << "\nGlobal memory (bytes): " << p.totalGlobalMem << "\nCUDA runtime: " << runtime
        << "\nExecuted kernel architecture: sm_" << architecture / 10;
    return out.str();
}
} // namespace astraflow

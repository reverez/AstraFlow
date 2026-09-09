#include "astraflow/core/cuda_euler1d.hpp"
#include "exact_riemann.hpp"
#include <catch2/catch_test_macros.hpp>
#include <iostream>
using namespace astraflow;
template <class T> void compare(double tolerance) {
    constexpr int n = 400;
    Euler1D<T> cpu(n);
    CudaEuler1D<T> gpu(n);
    cpu.sod();
    gpu.sod();
    cpu.advance(T(0.2));
    gpu.advance(T(0.2));
    auto a = cpu.state(), b = gpu.state();
    double maximum = 0, l1 = 0, l2 = 0;
    for (std::size_t k = 0; k < a.size(); ++k) {
        double d = double(a[k]) - double(b[k]);
        maximum = std::max(maximum, std::abs(d));
        l1 += std::abs(d) / a.size();
        l2 += d * d / a.size();
    }
    std::cout << "CUDA Sod FP" << sizeof(T) * 8 << " max=" << maximum << " L1=" << l1
              << " L2=" << std::sqrt(l2)
              << " CPU_ms/step=" << cpu.diagnostics().elapsed_ms / cpu.diagnostics().iterations
              << " GPU_ms/step=" << gpu.diagnostics().elapsed_ms / gpu.diagnostics().iterations
              << '\n';
    REQUIRE(maximum < tolerance);
    reference::Riemann exact({{1, 0, 0, 1}}, {{0.125, 0, 0, 0.1}}, 1.4);
    double error[4] = {};
    for (int i = 0; i < n; ++i) {
        State<T> u;
        for (int k = 0; k < 4; ++k)
            u[k] = b[k * n + i];
        auto w = primitive(u, Gas{});
        auto ref = exact.sample((i + 0.5) / n, 0.2);
        REQUIRE(physical(w, Gas{}));
        for (int k = 0; k < 4; ++k)
            error[k] += std::abs(double(w[k]) - ref[k]) / n;
    }
    REQUIRE(error[0] < 0.006);
    REQUIRE(error[1] < 0.008);
    REQUIRE(error[3] < 0.005);
    std::vector<State<T>> uniform(n, {{1, T(0.4), T(0.2), 1}});
    gpu.initialize(uniform);
    auto before = gpu.state();
    gpu.advance(T(0.1));
    auto after = gpu.state();
    for (std::size_t i = 0; i < before.size(); ++i)
        REQUIRE(std::abs(after[i] - before[i]) < tolerance);
}
TEST_CASE("CUDA Euler FP64 parity, Sod and uniform state", "[cuda]") { compare<double>(2e-11); }
TEST_CASE("CUDA Euler FP32 parity, Sod and uniform state", "[cuda]") { compare<float>(5e-5); }

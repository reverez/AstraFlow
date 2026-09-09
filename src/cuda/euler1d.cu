#include "astraflow/core/cuda_euler1d.hpp"
#include "memory/resources.cuh"
#include <cub/device/device_reduce.cuh>
#include <limits>
namespace astraflow {
namespace {
constexpr int block = 128;
template <class T>
__global__ void primitive_and_cfl(View<T> u, View<T> w, T *dt, Gas gas, T cfl, int *bad) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= u.size)
        return;
    auto p = primitive(u.get(i), gas);
    w.set(i, p);
    if (!physical(p, gas))
        atomicCAS(bad, 0, i + 1);
    dt[i] = cfl / (T(u.size) * (absolute(p[1]) + sound_speed(p, gas)));
}
template <class T>
__global__ void reconstruct_flux(View<T> w, View<T> f, Gas gas, Limiter limiter,
                                 unsigned long long *counts) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i > w.size)
        return;
    bool corrected = false;
    auto result = face_1d(w, i, gas, limiter, corrected);
    f.set(i, result.value);
    if (result.fallback)
        atomicAdd(counts, 1ULL);
    if (corrected)
        atomicAdd(counts + 1, 1ULL);
}
template <class T> __global__ void assemble_residual(View<T> f, View<T> r) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < r.size)
        r.set(i, T(r.size) * (f.get(i) - f.get(i + 1)));
}
template <class T>
__global__ void rk_update(View<T> base, View<T> from, View<T> out, View<T> r, T dt, bool second,
                          Gas gas, int *bad) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= out.size)
        return;
    auto next = from.get(i) + dt * r.get(i);
    if (second)
        next = T(0.5) * (base.get(i) + next);
    if (!physical(primitive(next, gas), gas))
        atomicCAS(bad, 0, i + 1);
    out.set(i, next);
}
} // namespace
template <class T> struct CudaEuler1D<T>::Impl {
    int n;
    Gas gas;
    T cfl;
    Limiter limiter;
    cuda::Buffer<T> u, w, stage, residual, flux, steps, minimum;
    cuda::Buffer<int> bad{1};
    cuda::Buffer<unsigned long long> counts{2};
    std::unique_ptr<cuda::Buffer<unsigned char>> scratch;
    std::size_t scratch_bytes = 0;
    Diagnostics diag;
    cuda::Event start, stop;
    Impl(int cells, Gas g, T c, Limiter l)
        : n(cells), gas(g), cfl(c), limiter(l), u(4 * cells), w(4 * cells), stage(4 * cells),
          residual(4 * cells), flux(4 * (cells + 1)), steps(cells), minimum(1) {
        AF_CUDA(cub::DeviceReduce::Min(nullptr, scratch_bytes, steps.data(), minimum.data(), n));
        scratch = std::make_unique<cuda::Buffer<unsigned char>>(scratch_bytes);
        bad.zero();
        counts.zero();
    }
    View<T> view(cuda::Buffer<T> &b, int size = 0) { return {b.data(), size ? size : n}; }
    void primitives(cuda::Buffer<T> &state) {
        primitive_and_cfl<<<(n + block - 1) / block, block>>>(view(state), view(w), steps.data(),
                                                              gas, cfl, bad.data());
        AF_CUDA(cudaGetLastError());
    }
    void rhs() {
        reconstruct_flux<<<(n + block) / block, block>>>(view(w), view(flux, n + 1), gas, limiter,
                                                         counts.data());
        AF_CUDA(cudaGetLastError());
        assemble_residual<<<(n + block - 1) / block, block>>>(view(flux, n + 1), view(residual));
        AF_CUDA(cudaGetLastError());
    }
    void check_state() {
        int failed = 0;
        AF_CUDA(cudaMemcpy(&failed, bad.data(), sizeof(int), cudaMemcpyDeviceToHost));
        if (failed)
            throw std::runtime_error("Invalid CUDA Euler state at cell " +
                                     std::to_string(failed - 1) + ", iteration " +
                                     std::to_string(diag.iterations));
    }
};
template <class T> CudaEuler1D<T>::CudaEuler1D(int n, Gas g, T c, Limiter l) {
    if (n < 4 || !(c > 0 && c <= T(0.5)) || !(g.gamma > 1) || !(g.gas_constant > 0))
        throw std::invalid_argument("Invalid CUDA Euler configuration");
    impl_ = std::make_unique<Impl>(n, g, c, l);
}
template <class T> CudaEuler1D<T>::~CudaEuler1D() = default;
template <class T> void CudaEuler1D<T>::initialize(const std::vector<State<T>> &w) {
    auto &s = *impl_;
    if (w.size() != std::size_t(s.n))
        throw std::invalid_argument("CUDA initial state size mismatch");
    std::vector<T> packed(4 * s.n);
    View<T> v{packed.data(), s.n};
    for (int i = 0; i < s.n; ++i) {
        if (!physical(w[i], s.gas))
            throw std::invalid_argument("Nonphysical CUDA initial state");
        v.set(i, conservative(w[i], s.gas));
    }
    s.u.upload(packed.data());
    s.bad.zero();
    s.counts.zero();
    s.diag = {};
}
template <class T> void CudaEuler1D<T>::sod() {
    std::vector<State<T>> w(impl_->n);
    for (int i = 0; i < impl_->n; ++i)
        w[i] = i < impl_->n / 2 ? State<T>{{1, 0, 0, 1}} : State<T>{{T(0.125), 0, 0, T(0.1)}};
    initialize(w);
}
template <class T> void CudaEuler1D<T>::advance(T end_time) {
    auto &s = *impl_;
    if (!std::isfinite(end_time) || end_time < T(s.diag.time))
        throw std::invalid_argument("Invalid CUDA end time");
    s.start.record();
    while (s.diag.time < double(end_time)) {
        s.primitives(s.u);
        AF_CUDA(cub::DeviceReduce::Min(s.scratch->data(), s.scratch_bytes, s.steps.data(),
                                       s.minimum.data(), s.n));
        T dt = 0;
        AF_CUDA(cudaMemcpy(&dt, s.minimum.data(), sizeof(T), cudaMemcpyDeviceToHost));
        s.check_state();
        dt = smaller(dt, T(double(end_time) - s.diag.time));
        if (!(dt > 0) || !std::isfinite(dt))
            throw std::runtime_error("Invalid CUDA timestep");
        s.rhs();
        rk_update<<<(s.n + block - 1) / block, block>>>(s.view(s.u), s.view(s.u), s.view(s.stage),
                                                        s.view(s.residual), dt, false, s.gas,
                                                        s.bad.data());
        AF_CUDA(cudaGetLastError());
        s.primitives(s.stage);
        s.rhs();
        rk_update<<<(s.n + block - 1) / block, block>>>(s.view(s.u), s.view(s.stage), s.view(s.u),
                                                        s.view(s.residual), dt, true, s.gas,
                                                        s.bad.data());
        AF_CUDA(cudaGetLastError());
        s.check_state();
        s.diag.time += double(dt);
        ++s.diag.iterations;
    }
    s.stop.record();
    s.diag.elapsed_ms += s.stop.since(s.start);
    auto counts = s.counts.download();
    s.diag.flux_fallbacks = counts[0];
    s.diag.reconstruction_fallbacks = counts[1];
}
template <class T> std::vector<T> CudaEuler1D<T>::state() const { return impl_->u.download(); }
template <class T> const Diagnostics &CudaEuler1D<T>::diagnostics() const { return impl_->diag; }
template class CudaEuler1D<float>;
template class CudaEuler1D<double>;
} // namespace astraflow

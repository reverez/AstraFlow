#include "astraflow/numerics/finite_volume.hpp"
#include "memory/resources.cuh"
#include <cub/device/device_reduce.cuh>
#include <limits>
namespace astraflow {
namespace {
constexpr int block_size = 128;
template <class T>
__global__ void fv_primitives(View<T> u, View<T> w, T *dt, const Cell *cells, const Face *faces,
                              Settings s, int *bad) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= u.size)
        return;
    auto p = primitive(u.get(i), s.gas);
    w.set(i, p);
    if (!physical(p, s.gas))
        atomicCAS(bad, 0, i + 1);
    dt[i] = fv_timestep(p, cells[i], faces, s);
}
template <class T>
__global__ void fv_gradients(View<T> w, View<T> sx, View<T> sr, const Cell *c, const Face *f,
                             Settings s) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < w.size)
        slopes(i, w, sx, sr, c, f, s);
}
template <class T>
__global__ void fv_faces(View<T> w, View<T> sx, View<T> sr, View<T> flux, const Face *f, Settings s,
                         unsigned long long *counts) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= flux.size)
        return;
    bool corrected = false;
    auto q = fv_flux(i, w, sx, sr, f, s, corrected);
    flux.set(i, q.value);
    if (q.fallback)
        atomicAdd(counts, 1ULL);
    if (corrected)
        atomicAdd(counts + 1, 1ULL);
}
template <class T>
__global__ void fv_assemble(View<T> f, View<T> res, View<T> w, const Cell *c, const Face *faces) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < res.size)
        res.set(i, fv_residual(i, f, w, c, faces));
}
template <class T>
__global__ void fv_rk(View<T> base, View<T> from, View<T> out, View<T> res, T dt, bool second,
                      Gas gas, int *bad, T *squares) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= out.size)
        return;
    auto r = res.get(i), next = from.get(i) + dt * r;
    if (second)
        next = T(0.5) * (base.get(i) + next);
    out.set(i, next);
    if (!physical(primitive(next, gas), gas))
        atomicCAS(bad, 0, i + 1);
    for (int k = 0; k < 4; ++k)
        squares[k * out.size + i] = r[k] * r[k];
}
template <class T> class GpuSolver final : public Solver {
    Settings s_;
    int n_, nf_, blocks_;
    bool ready_ = false;
    StepStats stats_;
    cuda::Buffer<Cell> cells_;
    cuda::Buffer<Face> faces_;
    cuda::Buffer<T> u_, w_, sx_, sr_, stage_, res_, flux_, dt_, minimum_, squares_, sums_;
    cuda::Buffer<int> bad_{1};
    cuda::Buffer<unsigned long long> counts_{2};
    std::unique_ptr<cuda::Buffer<unsigned char>> scratch_;
    std::size_t bytes_ = 0;
    cuda::Event start_, stop_;
    View<T> view(cuda::Buffer<T> &a, int size = 0) { return {a.data(), size ? size : n_}; }
    void check_state() {
        int bad = 0;
        AF_CUDA(cudaMemcpy(&bad, bad_.data(), sizeof(int), cudaMemcpyDeviceToHost));
        if (bad)
            throw std::runtime_error("Invalid CUDA FV state at cell " + std::to_string(bad - 1) +
                                     ", iteration " + std::to_string(stats_.iterations));
    }
    void convert(cuda::Buffer<T> &u) {
        fv_primitives<<<blocks_, block_size>>>(view(u), view(w_), dt_.data(), cells_.data(),
                                               faces_.data(), s_, bad_.data());
        AF_CUDA(cudaGetLastError());
    }
    void rhs() {
        fv_gradients<<<blocks_, block_size>>>(view(w_), view(sx_), view(sr_), cells_.data(),
                                              faces_.data(), s_);
        AF_CUDA(cudaGetLastError());
        fv_faces<<<(nf_ + block_size - 1) / block_size, block_size>>>(
            view(w_), view(sx_), view(sr_), view(flux_, nf_), faces_.data(), s_, counts_.data());
        AF_CUDA(cudaGetLastError());
        fv_assemble<<<blocks_, block_size>>>(view(flux_, nf_), view(res_), view(w_), cells_.data(),
                                             faces_.data());
        AF_CUDA(cudaGetLastError());
    }

  public:
    GpuSolver(const Mesh &m, Settings s)
        : s_(s), n_(int(m.cells.size())), nf_(int(m.faces.size())),
          blocks_((n_ + block_size - 1) / block_size), cells_(n_), faces_(nf_), u_(4 * n_),
          w_(4 * n_), sx_(4 * n_), sr_(4 * n_), stage_(4 * n_), res_(4 * n_), flux_(4 * nf_),
          dt_(n_), minimum_(1), squares_(4 * n_), sums_(4) {
        cells_.upload(m.cells.data());
        faces_.upload(m.faces.data());
        std::size_t min_bytes = 0, sum_bytes = 0;
        AF_CUDA(cub::DeviceReduce::Min(nullptr, min_bytes, dt_.data(), minimum_.data(), n_));
        AF_CUDA(cub::DeviceReduce::Sum(nullptr, sum_bytes, squares_.data(), sums_.data(), n_));
        bytes_ = larger(min_bytes, sum_bytes);
        scratch_ = std::make_unique<cuda::Buffer<unsigned char>>(bytes_);
        bad_.zero();
        counts_.zero();
    }
    void initialize(const std::vector<State<double>> &w) override {
        if (w.size() != std::size_t(n_))
            throw std::invalid_argument("CUDA FV initial size mismatch");
        std::vector<T> v(4 * n_);
        View<T> u{v.data(), n_};
        for (int i = 0; i < n_; ++i) {
            State<T> p;
            for (int k = 0; k < 4; ++k)
                p[k] = T(w[i][k]);
            if (!physical(p, s_.gas))
                throw std::invalid_argument("Nonphysical CUDA FV initialization");
            u.set(i, conservative(p, s_.gas));
        }
        u_.upload(v.data());
        bad_.zero();
        counts_.zero();
        stats_ = {};
        stats_.device_bytes = cells_.bytes() + faces_.bytes() +
                              (29 * n_ + 4 * nf_ + 5) * sizeof(T) + bad_.bytes() + counts_.bytes() +
                              bytes_;
        ready_ = true;
    }
    void step(double maximum_dt) override {
        if (!ready_ || !(maximum_dt > 0) || !std::isfinite(maximum_dt))
            throw std::invalid_argument("CUDA FV not ready or invalid time cap");
        start_.record();
        convert(u_);
        AF_CUDA(cub::DeviceReduce::Min(scratch_->data(), bytes_, dt_.data(), minimum_.data(), n_));
        T dt = 0;
        AF_CUDA(cudaMemcpy(&dt, minimum_.data(), sizeof(T), cudaMemcpyDeviceToHost));
        check_state();
        dt = smaller(dt, T(smaller(maximum_dt, double(std::numeric_limits<T>::max()))));
        if (!(dt > 0) || !std::isfinite(dt))
            throw std::runtime_error("Invalid CUDA FV timestep");
        rhs();
        fv_rk<<<blocks_, block_size>>>(view(u_), view(u_), view(stage_), view(res_), dt, false,
                                       s_.gas, bad_.data(), squares_.data());
        AF_CUDA(cudaGetLastError());
        convert(stage_);
        rhs();
        fv_rk<<<blocks_, block_size>>>(view(u_), view(stage_), view(u_), view(res_), dt, true,
                                       s_.gas, bad_.data(), squares_.data());
        AF_CUDA(cudaGetLastError());
        check_state();
        for (int k = 0; k < 4; ++k)
            AF_CUDA(cub::DeviceReduce::Sum(scratch_->data(), bytes_, squares_.data() + k * n_,
                                           sums_.data() + k, n_));
        auto sums = sums_.download();
        for (int k = 0; k < 4; ++k) {
            stats_.residual[k] = std::sqrt(double(sums[k]) / n_);
            if (!std::isfinite(stats_.residual[k]))
                throw std::runtime_error("Nonfinite CUDA FV residual");
        }
        auto counts = counts_.download();
        stats_.flux_fallbacks = counts[0];
        stats_.reconstruction_fallbacks = counts[1];
        stop_.record();
        stats_.iteration_ms = stop_.since(start_);
        stats_.total_ms += stats_.iteration_ms;
        stats_.dt = dt;
        stats_.time += dt;
        ++stats_.iterations;
    }
    std::vector<double> state() const override {
        auto v = u_.download();
        return {v.begin(), v.end()};
    }
    const StepStats &stats() const override { return stats_; }
};
} // namespace
std::unique_ptr<Solver> make_cuda_solver(const Mesh &m, Settings s, bool fp64) {
    m.validate();
    validate_settings(s);
    if (fp64)
        return std::make_unique<GpuSolver<double>>(m, s);
    return std::make_unique<GpuSolver<float>>(m, s);
}
} // namespace astraflow

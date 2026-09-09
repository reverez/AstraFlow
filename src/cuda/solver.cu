#include "astraflow/numerics/finite_volume.hpp"
#include "memory/resources.cuh"
#include <cub/device/device_reduce.cuh>
#include <cub/iterator/counting_input_iterator.cuh>
#include <cub/iterator/transform_input_iterator.cuh>
#include <limits>
namespace astraflow {
namespace {
constexpr int block_size = 128;
template <class T> struct SquaredResidual {
    View<T> residual;
    AF_HD State<T> operator()(int i) const {
        auto r = residual.get(i);
        for (int k = 0; k < 4; ++k)
            r[k] *= r[k];
        return r;
    }
};
template <class T> struct AddResidual {
    AF_HD State<T> operator()(State<T> a, State<T> b) const { return a + b; }
};
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
__global__ void fv_gradients(View<T> w, View<T> sx, View<T> sr, View<T> gx, View<T> gr,
                             const Cell *c, const Face *f, Settings s) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < w.size) {
        slopes(i, w, sx, sr, c, f, s);
        if (s.gas.viscosity > 0)
            transport_gradients(i, w, gx, gr, c, f, s);
    }
}
template <class T>
__global__ void fv_faces(View<T> w, View<T> sx, View<T> sr, View<T> gx, View<T> gr, View<T> flux,
                         const Cell *c, const Face *f, Settings s, unsigned long long *counts) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= flux.size)
        return;
    bool corrected = false;
    auto q = fv_flux(i, w, sx, sr, gx, gr, c, f, s, corrected);
    flux.set(i, q.value);
    if (q.fallback)
        atomicAdd(counts, 1ULL);
    if (corrected)
        atomicAdd(counts + 1, 1ULL);
}
template <class T>
__global__ void fv_assemble(View<T> f, View<T> res, View<T> w, View<T> gx, View<T> gr,
                            const Cell *c, const Face *faces, Settings s) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < res.size)
        res.set(i, fv_residual(i, f, w, gx, gr, c, faces, s));
}
template <class T>
__global__ void fv_rk(View<T> base, View<T> from, View<T> out, View<T> res, T dt, bool second,
                      Gas gas, int *bad) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= out.size)
        return;
    auto r = res.get(i), next = from.get(i) + dt * r;
    if (second)
        next = T(0.5) * (base.get(i) + next);
    out.set(i, next);
    if (!physical(primitive(next, gas), gas))
        atomicCAS(bad, 0, i + 1);
}
template <class T> class GpuSolver final : public Solver {
    Settings s_;
    int n_, nf_, blocks_;
    bool ready_ = false;
    StepStats stats_;
    cuda::Buffer<Cell> cells_;
    cuda::Buffer<Face> faces_;
    cuda::Buffer<T> u_, w_, sx_, sr_, gx_, gr_, stage_, res_, flux_, dt_, minimum_;
    cuda::Buffer<State<T>> sums_{1};
    cuda::Buffer<int> bad_{1};
    cuda::Buffer<unsigned long long> counts_{2};
    std::unique_ptr<cuda::Buffer<unsigned char>> scratch_;
    std::size_t bytes_ = 0;
    cuda::Event start_, stop_;
    std::array<cuda::Event, 24> profile_events_;
    int rk_stage_ = 0;
    void mark(int index) {
        if (s_.profile)
            profile_events_[index].record();
    }
    double interval(int begin, int end) {
        return profile_events_[end].since(profile_events_[begin]);
    }
    View<T> view(cuda::Buffer<T> &a, int size = 0) { return {a.data(), size ? size : n_}; }
    void reduce_residual(void *scratch, std::size_t &bytes) {
        using Iterator = cub::TransformInputIterator<State<T>, SquaredResidual<T>,
                                                     cub::CountingInputIterator<int>>;
        Iterator input(cub::CountingInputIterator<int>(0), SquaredResidual<T>{view(res_)});
        AF_CUDA(cub::DeviceReduce::Reduce(scratch, bytes, input, sums_.data(), n_, AddResidual<T>{},
                                          State<T>{}));
    }
    void check_state() {
        int bad = 0;
        AF_CUDA(cudaMemcpy(&bad, bad_.data(), sizeof(int), cudaMemcpyDeviceToHost));
        if (bad)
            throw std::runtime_error("Invalid CUDA FV state at cell " + std::to_string(bad - 1) +
                                     ", iteration " + std::to_string(stats_.iterations));
    }
    void convert(cuda::Buffer<T> &u) {
        mark(rk_stage_ ? 12 : 0);
        fv_primitives<<<blocks_, block_size>>>(view(u), view(w_), dt_.data(), cells_.data(),
                                               faces_.data(), s_, bad_.data());
        AF_CUDA(cudaGetLastError());
        mark(rk_stage_ ? 13 : 1);
    }
    void rhs() {
        int event = rk_stage_ ? 14 : 4;
        mark(event);
        fv_gradients<<<blocks_, block_size>>>(view(w_), view(sx_), view(sr_), view(gx_), view(gr_),
                                              cells_.data(), faces_.data(), s_);
        AF_CUDA(cudaGetLastError());
        mark(event + 1);
        mark(event + 2);
        fv_faces<<<(nf_ + block_size - 1) / block_size, block_size>>>(
            view(w_), view(sx_), view(sr_), view(gx_), view(gr_), view(flux_, nf_), cells_.data(),
            faces_.data(), s_, counts_.data());
        AF_CUDA(cudaGetLastError());
        mark(event + 3);
        mark(event + 4);
        fv_assemble<<<blocks_, block_size>>>(view(flux_, nf_), view(res_), view(w_), view(gx_),
                                             view(gr_), cells_.data(), faces_.data(), s_);
        AF_CUDA(cudaGetLastError());
        mark(event + 5);
    }

  public:
    GpuSolver(const Mesh &m, Settings s)
        : s_(s), n_(int(m.cells.size())), nf_(int(m.faces.size())),
          blocks_((n_ + block_size - 1) / block_size), cells_(n_), faces_(nf_), u_(4 * n_),
          w_(4 * n_), sx_(4 * n_), sr_(4 * n_), gx_(4 * n_), gr_(4 * n_), stage_(4 * n_),
          res_(4 * n_), flux_(4 * nf_), dt_(n_), minimum_(1) {
        cells_.upload(m.cells.data());
        faces_.upload(m.faces.data());
        std::size_t min_bytes = 0, sum_bytes = 0;
        AF_CUDA(cub::DeviceReduce::Min(nullptr, min_bytes, dt_.data(), minimum_.data(), n_));
        reduce_residual(nullptr, sum_bytes);
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
                              (33 * n_ + 4 * nf_ + 5) * sizeof(T) + bad_.bytes() + counts_.bytes() +
                              bytes_;
        ready_ = true;
    }
    void step(double maximum_dt) override {
        if (!ready_ || !(maximum_dt > 0) || !std::isfinite(maximum_dt))
            throw std::invalid_argument("CUDA FV not ready or invalid time cap");
        start_.record();
        rk_stage_ = 0;
        convert(u_);
        mark(2);
        AF_CUDA(cub::DeviceReduce::Min(scratch_->data(), bytes_, dt_.data(), minimum_.data(), n_));
        T dt = 0;
        AF_CUDA(cudaMemcpy(&dt, minimum_.data(), sizeof(T), cudaMemcpyDeviceToHost));
        check_state();
        mark(3);
        dt = smaller(dt, T(smaller(maximum_dt, double(std::numeric_limits<T>::max()))));
        if (!(dt > 0) || !std::isfinite(dt))
            throw std::runtime_error("Invalid CUDA FV timestep");
        rhs();
        mark(rk_stage_ ? 20 : 10);
        fv_rk<<<blocks_, block_size>>>(view(u_), view(u_), view(stage_), view(res_), dt, false,
                                       s_.gas, bad_.data());
        AF_CUDA(cudaGetLastError());
        mark(11);
        rk_stage_ = 1;
        convert(stage_);
        rhs();
        mark(rk_stage_ ? 20 : 10);
        fv_rk<<<blocks_, block_size>>>(view(u_), view(stage_), view(u_), view(res_), dt, true,
                                       s_.gas, bad_.data());
        AF_CUDA(cudaGetLastError());
        mark(21);
        mark(22);
        check_state();
        reduce_residual(scratch_->data(), bytes_);
        auto sums = sums_.download()[0];
        for (int k = 0; k < 4; ++k) {
            stats_.residual[k] = std::sqrt(double(sums[k]) / n_);
            if (!std::isfinite(stats_.residual[k]))
                throw std::runtime_error("Nonfinite CUDA FV residual");
        }
        auto counts = counts_.download();
        stats_.flux_fallbacks = counts[0];
        stats_.reconstruction_fallbacks = counts[1];
        mark(23);
        stop_.record();
        stats_.iteration_ms = stop_.since(start_);
        if (s_.profile)
            stats_.stage_ms = {interval(0, 1) + interval(12, 13),
                               interval(2, 3),
                               interval(4, 5) + interval(14, 15),
                               interval(6, 7) + interval(16, 17),
                               interval(8, 9) + interval(18, 19),
                               interval(10, 11) + interval(20, 21),
                               interval(22, 23)};
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

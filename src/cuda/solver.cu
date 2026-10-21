#include "astraflow/numerics/finite_volume.hpp"
#include "memory/resources.cuh"
#include <cstddef>
#include <cub/device/device_reduce.cuh>
#include <cub/iterator/counting_input_iterator.cuh>
#include <cub/iterator/transform_input_iterator.cuh>
#include <limits>
namespace astraflow {
namespace {
constexpr int block_size = 128, boundary_block = 256;
// Every per-step scalar diagnostic lives in one device record so a step needs a single
// device-to-host transfer. Fallback counters accumulate across the run; the rest is rewritten.
template <class T> struct alignas(16) Diagnostics {
    State<T> sums;
    T dt, dt_max, mass[2];
    int bad;
    unsigned long long counts[2];
};
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
// The global timestep stays on the device: the RK kernel reads the reduced minimum directly.
template <class T>
__global__ void fv_rk(View<T> base, View<T> from, View<T> out, View<T> res, const T *dt_device,
                      T cap, bool second, Gas gas, int *bad) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= out.size)
        return;
    T dt = smaller(*dt_device, cap);
    auto r = res.get(i), next = from.get(i) + dt * r;
    if (second)
        next = T(0.5) * (base.get(i) + next);
    out.set(i, next);
    if (!physical(primitive(next, gas), gas))
        atomicCAS(bad, 0, i + 1);
}
template <class T>
__global__ void pseudo_dt(View<T> w, T *local, const Cell *cells, const Face *faces, Settings s,
                          T cap, int *bad) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= w.size)
        return;
    T dt = smaller(cap, fv_timestep(w.get(i), cells[i], faces, s));
    local[i] = dt;
    if (!(dt > 0) || !std::isfinite(dt))
        atomicCAS(bad, 0, i + 1);
}
template <class T>
__global__ void pseudo_rk(View<T> base, View<T> from, View<T> out, View<T> res, const T *local,
                          const T *minimum, bool local_time, bool second, Gas gas, int *bad) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= out.size)
        return;
    T dt = local_time ? local[i] : *minimum;
    auto next = from.get(i) + dt * res.get(i);
    if (second)
        next = T(.5) * (base.get(i) + next);
    out.set(i, next);
    if (!physical(primitive(next, gas), gas))
        atomicCAS(bad, 0, i + 1);
}
// Deterministic single-block tree reduction of weighted boundary mass fluxes.
template <class T>
__global__ void boundary_mass(View<T> flux, const int *faces, const T *weights, int inlet,
                              int total, T *out) {
    __shared__ T partial[2][boundary_block];
    T in = 0, outflow = 0;
    for (int k = threadIdx.x; k < total; k += blockDim.x) {
        T v = weights[k] * flux.data[faces[k]];
        if (k < inlet)
            in += v;
        else
            outflow += v;
    }
    partial[0][threadIdx.x] = in;
    partial[1][threadIdx.x] = outflow;
    __syncthreads();
    for (int stride = blockDim.x / 2; stride > 0; stride /= 2) {
        if (threadIdx.x < stride) {
            partial[0][threadIdx.x] += partial[0][threadIdx.x + stride];
            partial[1][threadIdx.x] += partial[1][threadIdx.x + stride];
        }
        __syncthreads();
    }
    if (threadIdx.x == 0) {
        out[0] = partial[0][0];
        out[1] = partial[1][0];
    }
}
template <class T> class GpuSolver final : public Solver {
    Settings s_;
    SteadySettings steady_;
    std::unique_ptr<cuda::Buffer<T>> local_dt_;
    int n_, nf_, blocks_, inlet_faces_ = 0, boundary_faces_ = 0;
    bool ready_ = false;
    StepStats stats_;
    cuda::Buffer<Cell> cells_;
    cuda::Buffer<Face> faces_;
    cuda::Buffer<T> u_, w_, sx_, sr_, gx_, gr_, stage_, res_, flux_, dt_;
    cuda::Buffer<Diagnostics<T>> diag_{1};
    std::unique_ptr<cuda::Buffer<int>> boundary_index_;
    std::unique_ptr<cuda::Buffer<T>> boundary_weight_;
    std::unique_ptr<cuda::Buffer<unsigned char>> scratch_;
    std::size_t bytes_ = 0;
    cuda::Event start_, stop_;
    std::array<cuda::Event, 24> profile_events_;
    int rk_stage_ = 0;
    template <class F> F *field(std::size_t offset) {
        return reinterpret_cast<F *>(reinterpret_cast<unsigned char *>(diag_.data()) + offset);
    }
    State<T> *sums() { return field<State<T>>(offsetof(Diagnostics<T>, sums)); }
    T *dt_field() { return field<T>(offsetof(Diagnostics<T>, dt)); }
    T *dt_max_field() { return field<T>(offsetof(Diagnostics<T>, dt_max)); }
    T *mass_field() { return field<T>(offsetof(Diagnostics<T>, mass)); }
    int *bad_field() { return field<int>(offsetof(Diagnostics<T>, bad)); }
    unsigned long long *counts_field() {
        return field<unsigned long long>(offsetof(Diagnostics<T>, counts));
    }
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
        AF_CUDA(cub::DeviceReduce::Reduce(scratch, bytes, input, sums(), n_, AddResidual<T>{},
                                          State<T>{}));
    }
    void reduce_boundary() {
        if (boundary_faces_ == 0)
            return;
        boundary_mass<<<1, boundary_block>>>(view(flux_, nf_), boundary_index_->data(),
                                             boundary_weight_->data(), inlet_faces_,
                                             boundary_faces_, mass_field());
        AF_CUDA(cudaGetLastError());
    }
    // The only synchronising transfer of a step: residual sums, timestep, boundary fluxes,
    // invalid-state index and fallback counters arrive together.
    Diagnostics<T> fetch(const char *context) {
        auto d = diag_.download()[0];
        if (d.bad)
            throw std::runtime_error("Invalid CUDA FV state at cell " + std::to_string(d.bad - 1) +
                                     ", iteration " + std::to_string(stats_.iterations));
        for (int k = 0; k < 4; ++k) {
            stats_.residual[k] = std::sqrt(double(d.sums[k]) / n_);
            if (!std::isfinite(stats_.residual[k]))
                throw std::runtime_error(std::string("Nonfinite CUDA ") + context + " residual");
        }
        stats_.boundary_mass_flux = {double(d.mass[0]), double(d.mass[1])};
        stats_.flux_fallbacks = d.counts[0];
        stats_.reconstruction_fallbacks = d.counts[1];
        return d;
    }
    void convert(cuda::Buffer<T> &u) {
        mark(rk_stage_ ? 12 : 0);
        fv_primitives<<<blocks_, block_size>>>(view(u), view(w_), dt_.data(), cells_.data(),
                                               faces_.data(), s_, bad_field());
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
            faces_.data(), s_, counts_field());
        AF_CUDA(cudaGetLastError());
        mark(event + 3);
        mark(event + 4);
        fv_assemble<<<blocks_, block_size>>>(view(flux_, nf_), view(res_), view(w_), view(gx_),
                                             view(gr_), cells_.data(), faces_.data(), s_);
        AF_CUDA(cudaGetLastError());
        mark(event + 5);
    }
    T cap(double maximum_dt) const {
        return T(smaller(maximum_dt, double(std::numeric_limits<T>::max())));
    }
    void pseudo_step(double maximum_dt) {
        start_.record();
        rk_stage_ = 0;
        convert(u_);
        Settings timestep_settings = s_;
        timestep_settings.cfl = steady_.cfl;
        pseudo_dt<<<blocks_, block_size>>>(view(w_), local_dt_->data(), cells_.data(),
                                           faces_.data(), timestep_settings, cap(maximum_dt),
                                           bad_field());
        AF_CUDA(cudaGetLastError());
        AF_CUDA(
            cub::DeviceReduce::Min(scratch_->data(), bytes_, local_dt_->data(), dt_field(), n_));
        AF_CUDA(cub::DeviceReduce::Max(scratch_->data(), bytes_, local_dt_->data(), dt_max_field(),
                                       n_));
        rhs();
        pseudo_rk<<<blocks_, block_size>>>(view(u_), view(u_), view(stage_), view(res_),
                                           local_dt_->data(), dt_field(), steady_.local_time_step,
                                           false, s_.gas, bad_field());
        AF_CUDA(cudaGetLastError());
        convert(stage_);
        rhs();
        pseudo_rk<<<blocks_, block_size>>>(view(u_), view(stage_), view(u_), view(res_),
                                           local_dt_->data(), dt_field(), steady_.local_time_step,
                                           true, s_.gas, bad_field());
        AF_CUDA(cudaGetLastError());
        convert(u_);
        rhs();
        reduce_residual(scratch_->data(), bytes_);
        reduce_boundary();
        auto d = fetch("pseudo");
        stop_.record();
        stats_.iteration_ms = stop_.since(start_);
        stats_.total_ms += stats_.iteration_ms;
        stats_.pseudo_dt_min = d.dt;
        stats_.pseudo_dt_max = steady_.local_time_step ? d.dt_max : d.dt;
        ++stats_.iterations;
    }

  public:
    GpuSolver(const Mesh &m, Settings s, SteadySettings steady)
        : s_(s), steady_(steady), n_(int(m.cells.size())), nf_(int(m.faces.size())),
          blocks_((n_ + block_size - 1) / block_size), cells_(n_), faces_(nf_), u_(4 * n_),
          w_(4 * n_), sx_(4 * n_), sr_(4 * n_), gx_(4 * n_), gr_(4 * n_), stage_(4 * n_),
          res_(4 * n_), flux_(4 * nf_), dt_(n_) {
        cells_.upload(m.cells.data());
        faces_.upload(m.faces.data());
        std::vector<int> index;
        std::vector<T> weight;
        for (int pass = 0; pass < 2; ++pass)
            for (int f = 0; f < nf_; ++f) {
                const auto &face = m.faces[f];
                bool inlet = face.boundary == Boundary::Inlet,
                     outlet = face.boundary == Boundary::Outlet;
                if ((pass == 0 && inlet) || (pass == 1 && outlet)) {
                    index.push_back(f);
                    weight.push_back(T((inlet == (face.left < 0) ? 1 : -1) * face.area));
                }
                if (pass == 0 && f == nf_ - 1)
                    inlet_faces_ = int(index.size());
            }
        boundary_faces_ = int(index.size());
        if (boundary_faces_) {
            boundary_index_ = std::make_unique<cuda::Buffer<int>>(boundary_faces_);
            boundary_weight_ = std::make_unique<cuda::Buffer<T>>(boundary_faces_);
            boundary_index_->upload(index.data());
            boundary_weight_->upload(weight.data());
        }
        std::size_t min_bytes = 0, sum_bytes = 0;
        AF_CUDA(cub::DeviceReduce::Min(nullptr, min_bytes, dt_.data(), dt_field(), n_));
        reduce_residual(nullptr, sum_bytes);
        bytes_ = larger(min_bytes, sum_bytes);
        if (steady.enabled) {
            local_dt_ = std::make_unique<cuda::Buffer<T>>(n_);
            std::size_t max_bytes = 0;
            AF_CUDA(
                cub::DeviceReduce::Max(nullptr, max_bytes, local_dt_->data(), dt_max_field(), n_));
            bytes_ = larger(bytes_, max_bytes);
        }
        scratch_ = std::make_unique<cuda::Buffer<unsigned char>>(bytes_);
        diag_.zero();
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
        diag_.zero();
        stats_ = {};
        stats_.device_bytes = cells_.bytes() + faces_.bytes() + (33 * n_ + 4 * nf_) * sizeof(T) +
                              diag_.bytes() + bytes_;
        if (boundary_faces_)
            stats_.device_bytes += boundary_index_->bytes() + boundary_weight_->bytes();
        if (steady_.enabled)
            stats_.device_bytes += local_dt_->bytes();
        ready_ = true;
    }
    void step(double maximum_dt) override {
        if (!ready_ || !(maximum_dt > 0) || !std::isfinite(maximum_dt))
            throw std::invalid_argument("CUDA FV not ready or invalid time cap");
        if (steady_.enabled) {
            pseudo_step(maximum_dt);
            return;
        }
        start_.record();
        rk_stage_ = 0;
        convert(u_);
        mark(2);
        AF_CUDA(cub::DeviceReduce::Min(scratch_->data(), bytes_, dt_.data(), dt_field(), n_));
        mark(3);
        T limit = cap(maximum_dt);
        rhs();
        mark(10);
        fv_rk<<<blocks_, block_size>>>(view(u_), view(u_), view(stage_), view(res_), dt_field(),
                                       limit, false, s_.gas, bad_field());
        AF_CUDA(cudaGetLastError());
        mark(11);
        rk_stage_ = 1;
        convert(stage_);
        rhs();
        mark(20);
        fv_rk<<<blocks_, block_size>>>(view(u_), view(stage_), view(u_), view(res_), dt_field(),
                                       limit, true, s_.gas, bad_field());
        AF_CUDA(cudaGetLastError());
        mark(21);
        mark(22);
        reduce_residual(scratch_->data(), bytes_);
        reduce_boundary();
        auto d = fetch("FV");
        T dt = smaller(d.dt, limit);
        if (!(dt > 0) || !std::isfinite(dt))
            throw std::runtime_error("Invalid CUDA FV timestep");
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
std::unique_ptr<Solver> make_cuda_solver(const Mesh &m, Settings s, bool fp64,
                                         SteadySettings steady) {
    m.validate();
    validate_settings(s);
    steady.validate();
    if (fp64)
        return std::make_unique<GpuSolver<double>>(m, s, steady);
    return std::make_unique<GpuSolver<float>>(m, s, steady);
}
} // namespace astraflow

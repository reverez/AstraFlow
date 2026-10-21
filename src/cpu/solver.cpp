#include "astraflow/numerics/finite_volume.hpp"
#include <chrono>
#include <limits>
#include <stdexcept>
#ifdef _OPENMP
#include <omp.h>
#endif
namespace astraflow {
namespace {
int requested_threads = 1;
}
int cpu_threads() { return requested_threads; }
int cpu_hardware_threads() {
#ifdef _OPENMP
    return omp_get_num_procs();
#else
    return 1;
#endif
}
void set_cpu_threads(int n) {
    if (n < 1)
        throw std::invalid_argument("CPU thread count must be positive");
#ifndef _OPENMP
    if (n != 1)
        throw std::invalid_argument("This build has no OpenMP; only one CPU thread is available");
#endif
    requested_threads = n;
}
namespace {
// Fork/join costs dominate tiny meshes; below this many cells every loop stays serial.
constexpr int parallel_threshold = 2048;
// Boundary faces whose numerical mass flux defines the inlet/outlet balance.
struct BoundaryFace {
    int face;
    double weight; // inward for inlet, outward for outlet: area times orientation sign.
};
template <class T> class CpuSolver final : public Solver {
    Mesh mesh_;
    Settings settings_;
    SteadySettings steady_;
    std::vector<T> local_dt_;
    int n_, nf_;
    bool ready_ = false;
    std::vector<T> u_, w_, sx_, sr_, gx_, gr_, stage_, res_, flux_;
    std::vector<BoundaryFace> inlet_, outlet_;
    int threads_ = cpu_threads();
    bool parallel_ = threads_ > 1 && n_ >= parallel_threshold;
    StepStats stats_;
    View<T> view(std::vector<T> &v, int n = 0) { return {v.data(), n ? n : n_}; }
    void convert(std::vector<T> &source) {
        auto u = view(source), w = view(w_);
        int bad = -1;
#pragma omp parallel for schedule(static) if (parallel_) num_threads(threads_)
        for (int i = 0; i < n_; ++i) {
            auto p = primitive(u.get(i), settings_.gas);
            if (!physical(p, settings_.gas)) {
#pragma omp critical
                bad = bad < 0 ? i : smaller(bad, i);
            }
            w.set(i, p);
        }
        if (bad >= 0)
            throw std::runtime_error("Invalid CPU FV state at cell " + std::to_string(bad) +
                                     ", iteration " + std::to_string(stats_.iterations));
    }
    // Every loop below is embarrassingly parallel per cell/face; fallback counters reduce.
    void rhs() {
        auto w = view(w_), sx = view(sx_), sr = view(sr_), flux = view(flux_, nf_),
             res = view(res_), gx = view(gx_), gr = view(gr_);
        const Cell *cells = mesh_.cells.data();
        const Face *faces = mesh_.faces.data();
        bool viscous = settings_.gas.viscosity > 0;
#pragma omp parallel for schedule(static) if (parallel_) num_threads(threads_)
        for (int i = 0; i < n_; ++i) {
            slopes(i, w, sx, sr, cells, faces, settings_);
            if (viscous)
                transport_gradients(i, w, gx, gr, cells, faces, settings_);
        }
        std::uint64_t fallbacks = 0, corrections = 0;
#pragma omp parallel for schedule(static) if (parallel_) num_threads(threads_)                     \
    reduction(+ : fallbacks, corrections)
        for (int i = 0; i < nf_; ++i) {
            bool corrected = false;
            auto f = fv_flux(i, w, sx, sr, gx, gr, cells, faces, settings_, corrected);
            flux.set(i, f.value);
            fallbacks += f.fallback;
            corrections += corrected;
        }
        stats_.flux_fallbacks += fallbacks;
        stats_.reconstruction_fallbacks += corrections;
#pragma omp parallel for schedule(static) if (parallel_) num_threads(threads_)
        for (int i = 0; i < n_; ++i)
            res.set(i, fv_residual(i, flux, w, gx, gr, cells, faces, settings_));
    }
    // Serial, fixed-order reductions keep the reference statistics independent of thread count.
    void diagnostics() {
        stats_.residual = {};
        for (int k = 0; k < 4 * n_; ++k)
            stats_.residual[k / n_] += double(res_[k]) * double(res_[k]) / n_;
        for (auto &r : stats_.residual) {
            r = std::sqrt(r);
            if (!std::isfinite(r))
                throw std::runtime_error("Nonfinite FV residual");
        }
        stats_.boundary_mass_flux = {};
        for (const auto &b : inlet_)
            stats_.boundary_mass_flux[0] += b.weight * double(flux_[b.face]);
        for (const auto &b : outlet_)
            stats_.boundary_mass_flux[1] += b.weight * double(flux_[b.face]);
        stats_.threads = threads_;
    }
    T cap(double maximum_dt) const {
        return T(smaller(maximum_dt, double(std::numeric_limits<T>::max())));
    }
    void pseudo_step(double maximum_dt) {
        auto start = std::chrono::steady_clock::now();
        convert(u_);
        auto w = view(w_);
        Settings timestep_settings = settings_;
        timestep_settings.cfl = steady_.cfl;
        T limit = cap(maximum_dt), minimum = limit, maximum = 0;
        bool invalid = false;
#pragma omp parallel for schedule(static) if (parallel_) num_threads(threads_)                     \
    reduction(min : minimum) reduction(max : maximum) reduction(| : invalid)
        for (int i = 0; i < n_; ++i) {
            auto dt = smaller(limit, fv_timestep(w.get(i), mesh_.cells[i], mesh_.faces.data(),
                                                 timestep_settings));
            invalid |= !(dt > 0) || !std::isfinite(dt);
            local_dt_[i] = dt;
            minimum = smaller(minimum, dt);
            maximum = larger(maximum, dt);
        }
        if (invalid)
            throw std::runtime_error("Invalid local pseudo timestep");
        if (!steady_.local_time_step) {
            std::fill(local_dt_.begin(), local_dt_.end(), minimum);
            maximum = minimum;
        }
        rhs();
#pragma omp parallel for schedule(static) if (parallel_) num_threads(threads_)
        for (int k = 0; k < 4 * n_; ++k)
            stage_[k] = u_[k] + local_dt_[k % n_] * res_[k];
        convert(stage_);
        rhs();
#pragma omp parallel for schedule(static) if (parallel_) num_threads(threads_)
        for (int k = 0; k < 4 * n_; ++k)
            u_[k] = T(.5) * (u_[k] + stage_[k] + local_dt_[k % n_] * res_[k]);
        convert(u_);
        // Certify R(U) itself at the completed pseudo iterate, never dt_i R(U).
        rhs();
        diagnostics();
        stats_.pseudo_dt_min = minimum;
        stats_.pseudo_dt_max = maximum;
        ++stats_.iterations;
        stats_.iteration_ms =
            std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start)
                .count();
        stats_.total_ms += stats_.iteration_ms;
    }

  public:
    CpuSolver(const Mesh &m, Settings s, SteadySettings steady)
        : mesh_(m), settings_(s), steady_(steady), local_dt_(steady.enabled ? m.cells.size() : 0),
          n_(int(m.cells.size())), nf_(int(m.faces.size())), u_(4 * n_), w_(4 * n_), sx_(4 * n_),
          sr_(4 * n_), gx_(4 * n_), gr_(4 * n_), stage_(4 * n_), res_(4 * n_), flux_(4 * nf_) {
        for (int f = 0; f < nf_; ++f) {
            const auto &face = mesh_.faces[f];
            if (face.boundary == Boundary::Inlet)
                inlet_.push_back({f, (face.left < 0 ? 1.0 : -1.0) * face.area});
            else if (face.boundary == Boundary::Outlet)
                outlet_.push_back({f, (face.left < 0 ? -1.0 : 1.0) * face.area});
        }
    }
    void initialize(const std::vector<State<double>> &w) override {
        if (w.size() != std::size_t(n_))
            throw std::invalid_argument("FV initial state size mismatch");
        auto u = view(u_);
        for (int i = 0; i < n_; ++i) {
            State<T> p;
            for (int k = 0; k < 4; ++k)
                p[k] = T(w[i][k]);
            if (!physical(p, settings_.gas))
                throw std::invalid_argument("Nonphysical FV initial state");
            u.set(i, conservative(p, settings_.gas));
        }
        stats_ = {};
        stats_.threads = threads_;
        ready_ = true;
    }
    void step(double maximum_dt) override {
        if (!ready_ || !(maximum_dt > 0) || !std::isfinite(maximum_dt))
            throw std::invalid_argument("FV solver not ready or invalid time cap");
        if (steady_.enabled) {
            pseudo_step(maximum_dt);
            return;
        }
        auto start = std::chrono::steady_clock::now();
        convert(u_);
        T dt = cap(maximum_dt);
        auto w = view(w_);
#pragma omp parallel for schedule(static) if (parallel_) num_threads(threads_) reduction(min : dt)
        for (int i = 0; i < n_; ++i)
            dt = smaller(dt, fv_timestep(w.get(i), mesh_.cells[i], mesh_.faces.data(), settings_));
        if (!(dt > 0) || !std::isfinite(dt))
            throw std::runtime_error("Invalid FV timestep");
        rhs();
#pragma omp parallel for schedule(static) if (parallel_) num_threads(threads_)
        for (int k = 0; k < 4 * n_; ++k)
            stage_[k] = u_[k] + dt * res_[k];
        convert(stage_);
        rhs();
#pragma omp parallel for schedule(static) if (parallel_) num_threads(threads_)
        for (int k = 0; k < 4 * n_; ++k)
            u_[k] = T(0.5) * (u_[k] + stage_[k] + dt * res_[k]);
        convert(u_);
        diagnostics();
        stats_.dt = dt;
        stats_.time += dt;
        ++stats_.iterations;
        stats_.iteration_ms =
            std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start)
                .count();
        stats_.total_ms += stats_.iteration_ms;
    }
    std::vector<double> state() const override { return {u_.begin(), u_.end()}; }
    const StepStats &stats() const override { return stats_; }
};
} // namespace
std::unique_ptr<Solver> make_cpu_solver(const Mesh &m, Settings s, bool fp64,
                                        SteadySettings steady) {
    m.validate();
    validate_settings(s);
    steady.validate();
    if (fp64)
        return std::make_unique<CpuSolver<double>>(m, s, steady);
    return std::make_unique<CpuSolver<float>>(m, s, steady);
}
} // namespace astraflow

#include "astraflow/numerics/finite_volume.hpp"
#include <chrono>
#include <limits>
#include <stdexcept>
namespace astraflow {
namespace {
template <class T> class CpuSolver final : public Solver {
    Mesh mesh_;
    Settings settings_;
    int n_, nf_;
    bool ready_ = false;
    std::vector<T> u_, w_, sx_, sr_, gx_, gr_, stage_, res_, flux_;
    StepStats stats_;
    View<T> view(std::vector<T> &v, int n = 0) { return {v.data(), n ? n : n_}; }
    void convert(std::vector<T> &source) {
        auto u = view(source), w = view(w_);
        for (int i = 0; i < n_; ++i) {
            auto p = primitive(u.get(i), settings_.gas);
            if (!physical(p, settings_.gas))
                throw std::runtime_error("Invalid CPU FV state at cell " + std::to_string(i) +
                                         ", iteration " + std::to_string(stats_.iterations));
            w.set(i, p);
        }
    }
    void rhs() {
        auto w = view(w_), sx = view(sx_), sr = view(sr_), flux = view(flux_, nf_),
             res = view(res_), gx = view(gx_), gr = view(gr_);
        for (int i = 0; i < n_; ++i)
            slopes(i, w, sx, sr, mesh_.cells.data(), mesh_.faces.data(), settings_);
        if (settings_.gas.viscosity > 0)
            for (int i = 0; i < n_; ++i)
                transport_gradients(i, w, gx, gr, mesh_.cells.data(), mesh_.faces.data(),
                                    settings_);
        for (int i = 0; i < nf_; ++i) {
            bool corrected = false;
            auto f = fv_flux(i, w, sx, sr, gx, gr, mesh_.cells.data(), mesh_.faces.data(),
                             settings_, corrected);
            flux.set(i, f.value);
            stats_.flux_fallbacks += f.fallback;
            stats_.reconstruction_fallbacks += corrected;
        }
        for (int i = 0; i < n_; ++i)
            res.set(i, fv_residual(i, flux, w, gx, gr, mesh_.cells.data(), mesh_.faces.data(),
                                   settings_));
    }

  public:
    CpuSolver(const Mesh &m, Settings s)
        : mesh_(m), settings_(s), n_(int(m.cells.size())), nf_(int(m.faces.size())), u_(4 * n_),
          w_(4 * n_), sx_(4 * n_), sr_(4 * n_), gx_(4 * n_), gr_(4 * n_), stage_(4 * n_),
          res_(4 * n_), flux_(4 * nf_) {}
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
        ready_ = true;
    }
    void step(double maximum_dt) override {
        if (!ready_ || !(maximum_dt > 0) || !std::isfinite(maximum_dt))
            throw std::invalid_argument("FV solver not ready or invalid time cap");
        auto start = std::chrono::steady_clock::now();
        convert(u_);
        T dt = T(smaller(maximum_dt, double(std::numeric_limits<T>::max())));
        auto w = view(w_);
        for (int i = 0; i < n_; ++i)
            dt = smaller(dt, fv_timestep(w.get(i), mesh_.cells[i], mesh_.faces.data(), settings_));
        if (!(dt > 0) || !std::isfinite(dt))
            throw std::runtime_error("Invalid FV timestep");
        rhs();
        for (int k = 0; k < 4 * n_; ++k)
            stage_[k] = u_[k] + dt * res_[k];
        convert(stage_);
        rhs();
        stats_.residual = {};
        for (int k = 0; k < 4 * n_; ++k) {
            u_[k] = T(0.5) * (u_[k] + stage_[k] + dt * res_[k]);
            stats_.residual[k / n_] += double(res_[k]) * double(res_[k]) / n_;
        }
        convert(u_);
        for (double &x : stats_.residual) {
            x = std::sqrt(x);
            if (!std::isfinite(x))
                throw std::runtime_error("Nonfinite FV residual");
        }
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
std::unique_ptr<Solver> make_cpu_solver(const Mesh &m, Settings s, bool fp64) {
    m.validate();
    validate_settings(s);
    if (fp64)
        return std::make_unique<CpuSolver<double>>(m, s);
    return std::make_unique<CpuSolver<float>>(m, s);
}
} // namespace astraflow

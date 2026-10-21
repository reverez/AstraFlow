#pragma once
#include "astraflow/core/euler1d.hpp"
#include "astraflow/geometry/mesh.hpp"
#include <array>
#include <memory>
#include <string>
namespace astraflow {
struct Settings {
    Gas gas;
    double cfl = 0.4;
    Limiter limiter = Limiter::MC;
    bool no_slip = false, profile = false;
    double p0 = 1, t0 = 1, back_pressure = 0.1;
    double wall_temperature = 0, upper_wall_speed = 0;
};
struct SteadySettings {
    bool enabled = false, local_time_step = true;
    double cfl = 0.4;
    void validate() const;
};
struct StepStats {
    int iterations = 0;
    double time = 0, dt = 0, iteration_ms = 0, total_ms = 0;
    // No common physical clock exists for local pseudo-time. These are normalized step bounds.
    double pseudo_dt_min = 0, pseudo_dt_max = 0;
    std::array<double, 4> residual{};
    // Numerical inlet (inward) and outlet (outward) mass fluxes from the same spatial
    // evaluation that produced `residual`; nondimensional, per radian when axisymmetric.
    std::array<double, 2> boundary_mass_flux{};
    std::array<double, 7> stage_ms{};
    int threads = 1;
    std::uint64_t flux_fallbacks = 0, reconstruction_fallbacks = 0;
    std::size_t device_bytes = 0;
};
class Solver {
  public:
    virtual ~Solver() = default;
    virtual void initialize(const std::vector<State<double>> &primitives) = 0;
    virtual void step(double maximum_dt = 1e100) = 0;
    virtual std::vector<double> state() const = 0;
    virtual const StepStats &stats() const = 0;
};
std::unique_ptr<Solver> make_cpu_solver(const Mesh &, Settings, bool fp64 = true,
                                        SteadySettings = {});
#ifdef ASTRAFLOW_HAS_CUDA
std::unique_ptr<Solver> make_cuda_solver(const Mesh &, Settings, bool fp64 = false,
                                         SteadySettings = {});
#endif
void validate_settings(const Settings &);
// CPU worker threads for the reference backend. The default is one thread, which keeps the
// historical single-thread reference semantics; set_cpu_threads opts into OpenMP workers (bound
// to the solver at construction). cpu_hardware_threads reports the available logical CPUs.
int cpu_threads();
int cpu_hardware_threads();
void set_cpu_threads(int);
inline double boundary_mass_mismatch(const std::array<double, 2> &f) {
    double scale = f[0] < 0 ? -f[0] : f[0], o = f[1] < 0 ? -f[1] : f[1];
    scale = scale > o ? scale : o;
    double d = f[1] - f[0];
    return (d < 0 ? -d : d) / (scale > 1e-30 ? scale : 1e-30);
}
} // namespace astraflow

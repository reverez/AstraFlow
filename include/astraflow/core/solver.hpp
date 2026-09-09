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
    bool no_slip = false;
    double p0 = 1, t0 = 1, back_pressure = 0.1;
    double wall_temperature = 0, upper_wall_speed = 0;
};
struct StepStats {
    int iterations = 0;
    double time = 0, dt = 0, iteration_ms = 0, total_ms = 0;
    std::array<double, 4> residual{};
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
std::unique_ptr<Solver> make_cpu_solver(const Mesh &, Settings, bool fp64 = true);
#ifdef ASTRAFLOW_HAS_CUDA
std::unique_ptr<Solver> make_cuda_solver(const Mesh &, Settings, bool fp64 = false);
#endif
void validate_settings(const Settings &);
} // namespace astraflow

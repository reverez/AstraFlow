#include "astraflow/core/simulation.hpp"
#include "astraflow/analysis/engineering.hpp"
#include <algorithm>
#include <stdexcept>
namespace astraflow {
Simulation::Simulation(Config config)
    : config_(std::move(config)), monitor_(config_.convergence, config_.residual_tolerance) {
    config_.validate();
    auto s = config_.settings;
    bool nozzle = config_.problem == "rocket_nozzle" || config_.problem == "isentropic_nozzle";
    physical_mesh_ =
        nozzle ? nozzle_mesh(config_.geometry)
               : rectangular_mesh(config_.geometry.nx, config_.geometry.nr, config_.length,
                                  config_.height, config_.problem == "viscous_channel", false,
                                  config_.problem == "viscous_channel");
    scales_.length = nozzle ? config_.geometry.throat_radius : config_.height;
    scales_.temperature = s.t0;
    scales_.pressure = s.p0;
    scales_.density = s.p0 / (s.gas.gas_constant * s.t0);
    scales_.velocity = std::sqrt(s.gas.gas_constant * s.t0);
    scales_.time = scales_.length / scales_.velocity;
    Mesh m = physical_mesh_;
    double L = scales_.length;
    for (auto &p : m.vertices) {
        p.x /= L;
        p.r /= L;
    }
    for (auto &f : m.faces) {
        f.x /= L;
        f.r /= L;
        f.area /= m.axisymmetric ? L * L : L;
    }
    for (auto &c : m.cells) {
        c.x /= L;
        c.r /= L;
        c.dx /= L;
        c.dr /= L;
        c.volume /= m.axisymmetric ? L * L * L : L * L;
        c.planar_area /= L * L;
        c.radial_source *= L;
    }
    std::vector<State<double>> w;
    if (nozzle)
        w = nozzle_initial_state(physical_mesh_, config_.geometry, s);
    else
        for (auto c : physical_mesh_.cells) {
            if (config_.problem == "sod")
                w.push_back(c.x < config_.length / 2 ? State<double>{{1, 0, 0, 1}}
                                                     : State<double>{{0.125, 0, 0, 0.1}});
            else
                w.push_back(
                    {{scales_.density, 0.8 * s.upper_wall_speed * c.r / config_.height, 0, s.p0}});
        }
    for (auto &q : w) {
        q[0] /= scales_.density;
        q[1] /= scales_.velocity;
        q[2] /= scales_.velocity;
        q[3] /= scales_.pressure;
    }
    s.gas.gas_constant = 1;
    s.gas.viscosity /= scales_.density * scales_.velocity * L;
    s.gas.rho_floor /= scales_.density;
    s.gas.p_floor /= scales_.pressure;
    s.gas.temperature_floor /= scales_.temperature;
    s.p0 = 1;
    s.t0 = 1;
    s.back_pressure /= scales_.pressure;
    s.wall_temperature /= scales_.temperature;
    s.upper_wall_speed /= scales_.velocity;
    if (config_.backend == "cpu")
        solver_ = make_cpu_solver(m, s, config_.precision == "double");
#ifdef ASTRAFLOW_HAS_CUDA
    else
        solver_ = make_cuda_solver(m, s, config_.precision == "double");
#else
    else
        throw std::runtime_error("CUDA requested but this build has CUDA disabled");
#endif
    solver_->initialize(w);
}
void Simulation::step() {
    if (finished())
        return;
    double cap = config_.end_time > 0 ? (config_.end_time - stats().time) / scales_.time : 1e100;
    try {
        solver_->step(cap);
        const auto &s = solver_->stats();
        monitor_.initial(s.residual);
        if (s.iterations == 1 || s.iterations % config_.convergence.sampling_interval == 0) {
            engineering_ = engineering(physical_mesh_, state(), config_.settings.gas,
                                       config_.settings.back_pressure);
            sampled_iteration_ = s.iterations;
            monitor_.sample(s.iterations, s.residual, engineering_);
        }
    } catch (const std::exception &e) {
        failure_ = e.what();
        throw;
    }
}
bool Simulation::finished() const { return termination() != "running"; }
std::string Simulation::termination() const {
    const auto &s = solver_->stats();
    if (!failure_.empty())
        return "invalid_state";
    if (stopped_)
        return "user_stop";
    if (monitor_.converged())
        return "steady_converged";
    if (s.iterations >= config_.max_iterations)
        return "iteration_limit";
    if (config_.end_time > 0 && s.time * scales_.time >= config_.end_time * (1 - 1e-13))
        return "time_limit";
    return "running";
}
nlohmann::json Simulation::convergence_report() const {
    auto report = monitor_.report(solver_->stats().residual);
    report["converged"] = termination() == "steady_converged";
    report["termination_reason"] = termination();
    if (!failure_.empty())
        report["failure"] = failure_;
    return report;
}
std::vector<double> Simulation::state() const {
    auto u = solver_->state();
    std::size_t n = physical_mesh_.cells.size();
    double factors[4] = {scales_.density, scales_.density * scales_.velocity,
                         scales_.density * scales_.velocity, scales_.pressure};
    for (int k = 0; k < 4; ++k)
        for (std::size_t i = 0; i < n; ++i)
            u[k * n + i] *= factors[k];
    return u;
}
StepStats Simulation::stats() const {
    auto s = solver_->stats();
    s.time *= scales_.time;
    s.dt *= scales_.time;
    return s;
}
} // namespace astraflow

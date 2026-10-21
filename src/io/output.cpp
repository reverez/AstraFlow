#include "astraflow/io/output.hpp"
#include "astraflow/analysis/wall.hpp"
#include "astraflow/core/device.hpp"
#include <chrono>
#include <iomanip>
#include <stdexcept>
namespace astraflow {
namespace {
void json_file(const std::filesystem::path &p, const nlohmann::json &j) {
    std::ofstream out(p);
    out.exceptions(std::ios::failbit | std::ios::badbit);
    out << j.dump(2) << '\n';
}
} // namespace
void write_vtk(const std::filesystem::path &path, const Mesh &m, const std::vector<double> &u,
               Gas gas) {
    int n = int(m.cells.size());
    if (!u.empty() && u.size() != std::size_t(4 * n))
        throw std::invalid_argument("VTK state size mismatch");
    std::ofstream out(path);
    out.exceptions(std::ios::failbit | std::ios::badbit);
    out << std::setprecision(17);
    out << "<?xml version=\"1.0\"?>\n<VTKFile type=\"StructuredGrid\" version=\"0.1\" "
           "byte_order=\"LittleEndian\"><StructuredGrid WholeExtent=\"0 "
        << m.nx << " 0 " << m.nr << " 0 0\"><Piece Extent=\"0 " << m.nx << " 0 " << m.nr
        << " 0 0\">\n<Points><DataArray type=\"Float64\" NumberOfComponents=\"3\" "
           "format=\"ascii\">\n";
    for (auto p : m.vertices)
        out << p.x << ' ' << p.r << " 0\n";
    out << "</DataArray></Points>\n<CellData>\n";
    if (!u.empty())
        for (const char *name :
             {"density", "pressure", "temperature", "velocity", "Mach", "total_energy"}) {
            std::string field = name;
            out << "<DataArray type=\"Float64\" Name=\"" << name << "\" NumberOfComponents=\""
                << (field == "velocity" ? 3 : 1) << "\" format=\"ascii\">\n";
            for (int i = 0; i < n; ++i) {
                State<double> q;
                for (int k = 0; k < 4; ++k)
                    q[k] = u[k * n + i];
                auto w = primitive(q, gas);
                if (field == "velocity")
                    out << w[1] << ' ' << w[2] << " 0\n";
                else {
                    double v = field == "density"       ? w[0]
                               : field == "pressure"    ? w[3]
                               : field == "temperature" ? w[3] / (w[0] * gas.gas_constant)
                               : field == "Mach" ? std::hypot(w[1], w[2]) / sound_speed(w, gas)
                                                 : q[3] / q[0];
                    out << v << '\n';
                }
            }
            out << "</DataArray>\n";
        }
    out << "</CellData></Piece></StructuredGrid></VTKFile>\n";
}
nlohmann::json export_state(const Simulation &sim) {
    auto state = sim.state();
    int n = int(sim.mesh().cells.size());
    nlohmann::json rows = nlohmann::json::array();
    for (int i = 0; i < n; ++i)
        rows.push_back({state[i], state[n + i], state[2 * n + i], state[3 * n + i]});
    auto effective = sim.config().json();
    nlohmann::json configuration;
    for (auto key : {"problem", "gas", "geometry", "mesh", "boundary_conditions", "numerics"})
        configuration[key] = effective.at(key);
    return {{"schema", 1},
            {"astraflow_version", version()},
            {"iterations", sim.stats().iterations},
            {"physical_simulated_time", sim.stats().time},
            {"time_integration", sim.time_integration()},
            {"layout", "row-major j*nx+i; [rho, rho*u, rho*v, rho*E] in configuration units"},
            {"configuration", configuration},
            {"physical_conservative", rows}};
}
std::vector<State<double>> import_state(const nlohmann::json &initialization,
                                        const Config &config) {
    if (initialization.at("schema") != 1)
        throw std::invalid_argument("Unknown initial-state schema");
    auto effective = config.json();
    for (auto key : {"problem", "gas", "geometry", "mesh", "boundary_conditions", "numerics"})
        if (initialization.at("configuration").at(key) != effective.at(key))
            throw std::invalid_argument("Initial-state configuration mismatch");
    std::vector<State<double>> initial;
    for (auto row : initialization.at("physical_conservative")) {
        if (!row.is_array() || row.size() != 4)
            throw std::invalid_argument("Invalid initial conservative row");
        State<double> q;
        for (int k = 0; k < 4; ++k)
            q[k] = row[k].get<double>();
        auto w = primitive(q, config.settings.gas);
        if (!physical(w, config.settings.gas))
            throw std::invalid_argument("Nonphysical initial state");
        initial.push_back(w);
    }
    if (initial.size() != std::size_t(config.geometry.nx * config.geometry.nr))
        throw std::invalid_argument("Initial-state grid size mismatch");
    return initial;
}
RunOutput::RunOutput(const Simulation &sim) {
    directory_ = sim.config().output;
    if (directory_.empty())
        directory_ =
            std::filesystem::path("runs") /
            ("run-" + std::to_string(std::chrono::system_clock::now().time_since_epoch().count()));
    if (std::filesystem::exists(directory_) && !std::filesystem::is_empty(directory_))
        throw std::runtime_error("Output directory is not empty; choose a new path: " +
                                 directory_.string());
    std::filesystem::create_directories(directory_);
    auto config = sim.config().json();
    config["output"]["directory"] = directory_.string();
    json_file(directory_ / "config.json", config);
    write_vtk(directory_ / "mesh.vts", sim.mesh(), {}, sim.config().settings.gas);
    csv_.open(directory_ / "convergence.csv");
    csv_.exceptions(std::ios::failbit | std::ios::badbit);
    csv_ << std::setprecision(17)
         << "iteration,time,dt,mass_residual,axial_momentum_residual,radial_momentum_residual,"
            "energy_residual,iteration_ms,pseudo_dt_min_normalized,pseudo_dt_max_normalized,wall_"
            "elapsed_ms\n";
    engineering_csv_.open(directory_ / "engineering.csv");
    engineering_csv_.exceptions(std::ios::failbit | std::ios::badbit);
    engineering_csv_ << std::setprecision(17) << "iteration,time";
    for (auto field : stability_fields)
        engineering_csv_ << ',' << field;
    engineering_csv_
        << ",mass_conservation_error,mass_flow_spread,flux_fallbacks,reconstruction_fallbacks,"
           "interior_inlet_exit_mass_estimate_mismatch,numerical_boundary_mass_flux_mismatch\n";
    balance_history_.open(directory_ / "balance.jsonl");
    balance_history_.exceptions(std::ios::failbit | std::ios::badbit);
    json_file(directory_ / "initial_balance.json", sim.balance());
    json_file(directory_ / "metadata.json",
              {{"time_integration", sim.time_integration()},
               {"physical_time_columns", {"time", "dt"}},
               {"pseudo_time_coordinate", "none: local timesteps have no common elapsed time"},
               {"pseudo_timestep_units", "nondimensional"},
               {"residual_sampling",
                sim.config().steady.enabled ? "completed_pseudo_iterate" : "second_ssprk2_stage"},
               {"convergence_mass_gate", sim.config().convergence.mass_gate_key()},
               {"numerical_boundary_mass_flux_sampling",
                sim.config().steady.enabled ? "completed_pseudo_iterate" : "second_ssprk2_stage"},
               {"balance_interval", sim.config().convergence.balance_interval},
               {"astraflow_version", version()}});
    started_ = std::chrono::steady_clock::now();
}
void RunOutput::record(const Simulation &sim) {
    auto s = sim.stats();
    if (sim.sampled_iteration() == s.iterations && s.iterations != last_engineering_iteration_ &&
        s.iterations > 0) {
        last_engineering_iteration_ = s.iterations;
        engineering_csv_ << s.iterations << ',' << s.time;
        const auto &e = sim.sampled_engineering();
        for (auto field : stability_fields)
            engineering_csv_ << ',' << e.at(field);
        engineering_csv_ << ',' << e.at("mass_conservation_error") << ','
                         << e.at("mass_flow_spread") << ',' << s.flux_fallbacks << ','
                         << s.reconstruction_fallbacks << ','
                         << e.at("interior_inlet_exit_mass_estimate_mismatch") << ','
                         << e.at("numerical_boundary_mass_flux_mismatch") << '\n';
        if (sim.sampled_balance_iteration() == s.iterations &&
            s.iterations != last_balance_iteration_) {
            last_balance_iteration_ = s.iterations;
            auto b = sim.sampled_balance();
            b["iteration"] = s.iterations;
            b["physical_simulated_time"] = s.time;
            b["time_integration"] = sim.time_integration();
            b["interior_inlet_exit_mass_estimate_mismatch"] =
                e.at("interior_inlet_exit_mass_estimate_mismatch");
            balance_history_ << b.dump() << '\n';
            balance_history_.flush();
        }
        engineering_csv_.flush();
    }
    if (s.iterations == last_iteration_)
        return;
    last_iteration_ = s.iterations;
    csv_ << s.iterations << ',' << s.time << ',' << s.dt;
    for (double r : s.residual)
        csv_ << ',' << r;
    csv_ << ',' << s.iteration_ms << ',' << s.pseudo_dt_min << ',' << s.pseudo_dt_max << ','
         << std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - started_)
                .count()
         << '\n';
    csv_.flush();
}
nlohmann::json RunOutput::finish(const Simulation &sim, double wall_ms) {
    record(sim);
    csv_.flush();
    engineering_csv_.flush();
    auto result = sim.sampled_engineering();
    if (sim.termination() != "invalid_state") {
        auto state = sim.state();
        write_vtk(directory_ / "final_state.vts", sim.mesh(), state, sim.config().settings.gas);
        result = engineering(sim.mesh(), state, sim.config().settings.gas,
                             sim.config().settings.back_pressure);
        auto b = sim.balance();
        result["numerical_boundary_mass_flux_mismatch"] =
            b.at("numerical_boundary_mass_flux_mismatch");
        json_file(directory_ / "final_balance.json", b);
        if (sim.config().settings.gas.viscosity > 0 && sim.config().settings.no_slip) {
            auto samples = wall_distribution(sim.mesh(), state, sim.config().settings);
            if (!samples.empty()) {
                result["wall"] = wall_summary(sim.mesh(), samples);
                std::ofstream wall(directory_ / "wall.csv");
                wall.exceptions(std::ios::failbit | std::ios::badbit);
                wall << std::setprecision(17)
                     << "x,r,face_length,face_area,traction_x,traction_r,heat_flux_into_wall,"
                        "pressure,temperature,upper_wall\n";
                for (const auto &w : samples)
                    wall << w.x << ',' << w.r << ',' << w.length << ',' << w.area << ','
                         << w.traction_x << ',' << w.traction_r << ',' << w.heat_flux_into_wall
                         << ',' << w.pressure << ',' << w.temperature << ',' << w.upper << '\n';
            }
        }
    }
    result.update(sim.convergence_report());
    auto s = sim.stats();
    result["backend"] = sim.config().backend;
    result["precision"] = sim.config().precision;
    result["grid"] = {sim.mesh().nx, sim.mesh().nr};
    result["iterations"] = s.iterations;
    result["simulated_time"] = s.time;
    result["physical_simulated_time"] = s.time;
    result["time_integration"] = sim.time_integration();
    result["pseudo_iterations"] = sim.config().steady.enabled ? s.iterations : 0;
    result["pseudo_dt_min_normalized"] = s.pseudo_dt_min;
    result["pseudo_dt_max_normalized"] = s.pseudo_dt_max;
    result["convergence_mass_gate"] = sim.config().convergence.mass_gate_key();
    result["astraflow_version"] = version();
    result["cpu_threads"] = s.threads;
    result["termination"] = sim.termination();
    result["residuals_nondimensional"] = s.residual;
    result["flux_fallbacks"] = s.flux_fallbacks;
    result["reconstruction_fallbacks"] = s.reconstruction_fallbacks;
    json_file(directory_ / "summary.json", result);
    const auto &ref = sim.scales();
    json_file(directory_ / "performance.json",
              {{"time_integration", sim.time_integration()},
               {"wall_ms", wall_ms},
               {"solver_ms", s.total_ms},
               {"mean_iteration_ms", s.iterations ? s.total_ms / s.iterations : 0},
               {"device_bytes", s.device_bytes},
               {"cells", sim.mesh().cells.size()},
               {"reference_scales",
                {{"length", ref.length},
                 {"density", ref.density},
                 {"velocity", ref.velocity},
                 {"temperature", ref.temperature},
                 {"pressure", ref.pressure},
                 {"time", ref.time}}}});
    return result;
}
} // namespace astraflow

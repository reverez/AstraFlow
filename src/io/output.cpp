#include "astraflow/io/output.hpp"
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
            "energy_residual,iteration_ms\n";
}
void RunOutput::record(const Simulation &sim) {
    auto s = sim.stats();
    if (s.iterations == last_iteration_)
        return;
    last_iteration_ = s.iterations;
    csv_ << s.iterations << ',' << s.time << ',' << s.dt;
    for (double r : s.residual)
        csv_ << ',' << r;
    csv_ << ',' << s.iteration_ms << '\n';
}
nlohmann::json RunOutput::finish(const Simulation &sim, double wall_ms) {
    record(sim);
    csv_.flush();
    auto state = sim.state();
    write_vtk(directory_ / "final_state.vts", sim.mesh(), state, sim.config().settings.gas);
    auto result = engineering(sim.mesh(), state, sim.config().settings.gas,
                              sim.config().settings.back_pressure);
    auto s = sim.stats();
    result["backend"] = sim.config().backend;
    result["precision"] = sim.config().precision;
    result["grid"] = {sim.mesh().nx, sim.mesh().nr};
    result["iterations"] = s.iterations;
    result["simulated_time"] = s.time;
    result["termination"] = sim.termination();
    result["residuals_nondimensional"] = s.residual;
    result["flux_fallbacks"] = s.flux_fallbacks;
    result["reconstruction_fallbacks"] = s.reconstruction_fallbacks;
    json_file(directory_ / "summary.json", result);
    const auto &ref = sim.scales();
    json_file(directory_ / "performance.json",
              {{"wall_ms", wall_ms},
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

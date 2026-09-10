#include "astraflow/analysis/balance.hpp"
#include "astraflow/analysis/engineering.hpp"
#include "astraflow/numerics/finite_volume.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <set>
#include <sstream>
using namespace astraflow;
namespace {
void write_json(const std::filesystem::path &p, const nlohmann::json &j) {
    if (std::filesystem::exists(p))
        throw std::runtime_error("Output exists: " + p.string());
    std::ofstream f(p);
    f.exceptions(std::ios::badbit | std::ios::failbit);
    f << j.dump(2) << '\n';
}
int branch(double a, double b, Limiter limiter) {
    if (limiter == Limiter::FirstOrderDiagnostic)
        return 5;
    if (a * b <= 0)
        return 0;
    if (limiter == Limiter::VanLeer)
        return 4;
    double c = .5 * (a + b);
    if (std::abs(c) <= 2 * std::min(std::abs(a), std::abs(b)))
        return 1;
    return std::abs(a) < std::abs(b) ? 2 : 3;
}
} // namespace
int main(int argc, char **argv) {
    try {
        if (argc == 6 && std::string(argv[1]) == "prolong") {
            auto c = Config::read(argv[2]), f = Config::read(argv[4]);
            auto cj = c.json(), fj = f.json();
            for (auto key : {"problem", "gas", "geometry", "boundary_conditions", "numerics"})
                if (cj[key] != fj[key])
                    throw std::invalid_argument("Prolongation physical configuration differs");
            std::ifstream in(argv[3]);
            nlohmann::json source;
            in >> source;
            auto rows = source.at("physical_conservative");
            int n = c.geometry.nx * c.geometry.nr;
            if (rows.size() != std::size_t(n))
                throw std::invalid_argument("Parent state size mismatch");
            std::vector<double> u(4 * n);
            for (int i = 0; i < n; ++i)
                for (int k = 0; k < 4; ++k)
                    u[k * n + i] = rows[i][k].get<double>();
            nlohmann::json audit;
            auto state = conservative_prolong(nozzle_mesh(c.geometry), nozzle_mesh(f.geometry), u,
                                              c.settings.gas, audit);
            nlohmann::json output = nlohmann::json::array();
            for (auto q : state)
                output.push_back({q[0], q[1], q[2], q[3]});
            write_json(argv[5], {{"schema", 1},
                                 {"configuration", fj},
                                 {"physical_conservative", output},
                                 {"prolongation", audit},
                                 {"source", source.value("provenance", nlohmann::json::object())}});
            std::cout << audit.dump(2) << '\n';
            return 0;
        }
        if (argc != 8)
            throw std::invalid_argument(
                "Usage: residual_probe normalized-config primitives duration CFL output-directory "
                "mc|van_leer|first_order sample-steps; or residual_probe prolong coarse-config "
                "conservative-json fine-config output-json");
        auto c = Config::read(argv[1]);
        auto mesh = nozzle_mesh(c.geometry);
        Settings s = c.settings;
        s.cfl = std::stod(argv[4]);
        std::string mode = argv[6];
        if (mode == "mc")
            s.limiter = Limiter::MC;
        else if (mode == "van_leer")
            s.limiter = Limiter::VanLeer;
        else if (mode == "first_order")
            s.limiter = Limiter::FirstOrderDiagnostic;
        else
            throw std::invalid_argument("Unknown diagnostic mode");
        double duration = std::stod(argv[3]);
        int interval = std::stoi(argv[7]);
        if (!(duration > 0) || !std::isfinite(duration) || interval < 1)
            throw std::invalid_argument("Invalid diagnostic duration/sampling");
        std::filesystem::path dir = argv[5];
        if (std::filesystem::exists(dir))
            throw std::runtime_error("Diagnostic directory exists");
        std::filesystem::create_directories(dir);
        std::vector<State<double>> initial(mesh.cells.size());
        std::ifstream input(argv[2]);
        for (auto &q : initial)
            for (int k = 0; k < 4; ++k)
                if (!(input >> q[k]))
                    throw std::invalid_argument("Invalid primitive input");
        auto solver = make_cuda_solver(mesh, s, true);
        solver->initialize(initial);
        Settings ref = s;
        ref.limiter = Limiter::MC;
        auto reference = finite_volume_balance(mesh, solver->state(), ref);
        std::set<int> watch;
        std::vector<int> sorted(mesh.cells.size());
        std::iota(sorted.begin(), sorted.end(), 0);
        for (int k = 0; k < 4; ++k) {
            std::sort(sorted.begin(), sorted.end(), [&](int a, int b) {
                return std::abs(reference.total[a][k]) > std::abs(reference.total[b][k]);
            });
            for (int j = 0; j < 8; ++j)
                watch.insert(sorted[j]);
        }
        for (double x : {.15, .43}) {
            int i = std::clamp(int(x * mesh.nx), 0, mesh.nx - 1);
            for (int j : {0, 1, 2, 3, mesh.nr / 2, mesh.nr - 1})
                watch.insert(j * mesh.nx + i);
        }
        write_json(
            dir / "metadata.json",
            {{"diagnostic_only", true},
             {"mode", mode},
             {"cfl", s.cfl},
             {"normalized_duration", duration},
             {"sampling_steps", interval},
             {"watch_cells", watch},
             {"configuration", c.json()},
             {"balance_convention",
              "integral R = -sum outward(Fconv-Fvisc)A + integral S; per radian, nondimensional"}});
        std::ofstream history(dir / "balance.jsonl"), tracks(dir / "tracked.csv");
        history.exceptions(std::ios::failbit | std::ios::badbit);
        tracks.exceptions(std::ios::failbit | std::ios::badbit);
        tracks << std::setprecision(17);
        std::ostringstream header;
        header << "iteration,time,cell,i,j,x_over_L,r,rho,u,v,p";
        for (std::string prefix :
             {"rhs", "convective", "transport", "source", "sx", "sr", "gx", "gr", "left_jump",
              "right_jump", "bottom_jump", "top_jump", "branch_x", "branch_r"}) {
            std::array<const char *, 4> fields = {"rho", "u", "v", "p"};
            if (prefix == "rhs" || prefix == "convective" || prefix == "transport" ||
                prefix == "source")
                fields = {"mass", "axial", "radial", "energy"};
            if (prefix == "gx" || prefix == "gr")
                fields[3] = "T";
            for (auto field : fields)
                header << ',' << prefix << '_' << field;
        }
        header << ",reconstruction_faces,fallback_faces\n";
        tracks << header.str();
        auto sample = [&](bool dump) {
            auto u = solver->state();
            auto b = finite_volume_balance(mesh, u, s);
            auto st = solver->stats();
            auto g = b.global;
            g["iteration"] = st.iterations;
            g["time"] = st.time;
            g["rk_stage_rms"] = st.residual;
            g["interior_cell_estimate"] = engineering(mesh, u, s.gas, s.back_pressure);
            g["flux_fallbacks"] = st.flux_fallbacks;
            g["reconstruction_fallbacks"] = st.reconstruction_fallbacks;
            history << g.dump() << '\n';
            history.flush();
            int n = int(mesh.cells.size());
            std::vector<double> w(4 * n);
            for (int i = 0; i < n; ++i)
                for (int k = 0; k < 4; ++k)
                    w[k * n + i] = b.primitive[i][k];
            View<double> wv{w.data(), n};
            auto emit = [&](std::ostream &out, int i) {
                auto cell = mesh.cells[i];
                out << st.iterations << ',' << st.time << ',' << i << ',' << i % mesh.nx << ','
                    << i / mesh.nx << ',' << cell.x / c.geometry.length() << ',' << cell.r;
                for (int k = 0; k < 4; ++k)
                    out << ',' << b.primitive[i][k];
                State<double> jumps[4];
                jumps[0] = b.primitive[i] - neighbor(wv, cell, 0, mesh.faces.data(), s);
                jumps[1] = neighbor(wv, cell, 1, mesh.faces.data(), s) - b.primitive[i];
                jumps[2] = b.primitive[i] - neighbor(wv, cell, 2, mesh.faces.data(), s);
                jumps[3] = neighbor(wv, cell, 3, mesh.faces.data(), s) - b.primitive[i];
                for (auto q : {b.total[i], b.convective[i], b.transport[i], b.source[i], b.sx[i],
                               b.sr[i], b.gx[i], b.gr[i], jumps[0], jumps[1], jumps[2], jumps[3]})
                    for (int k = 0; k < 4; ++k)
                        out << ',' << q[k];
                for (int d = 0; d < 2; ++d)
                    for (int k = 0; k < 4; ++k)
                        out << ',' << branch(jumps[2 * d][k], jumps[2 * d + 1][k], s.limiter);
                out << ',' << b.reconstruction_faces[i] << ',' << b.fallback_faces[i] << '\n';
            };
            for (int i : watch)
                emit(tracks, i);
            tracks.flush();
            if (dump) {
                std::ofstream out(dir / (st.iterations == 0 ? "initial.csv" : "final.csv"));
                out << std::setprecision(17) << header.str();
                for (int i = 0; i < n; ++i)
                    emit(out, i);
            }
            return g;
        };
        auto first = sample(true);
        while (solver->stats().time < duration) {
            solver->step(duration - solver->stats().time);
            if (solver->stats().iterations % interval == 0 && solver->stats().time < duration)
                sample(false);
        }
        auto last = sample(true);
        write_json(dir / "summary.json", {{"diagnostic_only", true},
                                          {"initial", first},
                                          {"final", last},
                                          {"iterations", solver->stats().iterations},
                                          {"normalized_duration", solver->stats().time}});
        std::cout << mode << " steps=" << solver->stats().iterations
                  << " time=" << solver->stats().time << " initial/final RMS=" << first["rms"]
                  << ' ' << last["rms"] << '\n';
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}

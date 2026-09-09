// Targeted steady-residual diagnostic. Calls the existing shared FV mathematics unchanged.
#include "astraflow/io/config.hpp"
#include "astraflow/numerics/finite_volume.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace astraflow;
int main(int argc, char **argv) {
    try {
        if (argc != 6)
            throw std::invalid_argument(
                "Usage: diagnose_closure normalized-config primitives.txt steps CFL output-prefix");
        auto config = Config::read(argv[1]);
        auto mesh = nozzle_mesh(config.geometry);
        auto settings = config.settings;
        settings.cfl = std::stod(argv[4]);
        int steps = std::stoi(argv[3]), n = int(mesh.cells.size()), nf = int(mesh.faces.size());
        std::ifstream input(argv[2]);
        std::vector<State<double>> initial(n);
        for (auto &w : initial)
            for (int k = 0; k < 4; ++k)
                if (!(input >> w[k]))
                    throw std::runtime_error("Invalid primitive input");
        std::string prefix = argv[5];
        for (const auto &suffix : {"-history.csv", "-initial.csv", "-final.csv"})
            if (std::filesystem::exists(prefix + suffix))
                throw std::runtime_error("Diagnostic output already exists; choose a new prefix");
        auto solver = make_cuda_solver(mesh, settings, true);
        solver->initialize(initial);
        std::ofstream log(prefix + "-history.csv");
        log << std::setprecision(17) << "iteration,time,mass,axial,radial,energy\n";
        auto diagnose = [&](const std::string &suffix) {
            auto u = solver->state();
            std::vector<double> w(4 * n), sx(4 * n), sr(4 * n), gx(4 * n), gr(4 * n), flux(4 * nf);
            View<double> uv{u.data(), n}, wv{w.data(), n}, sxv{sx.data(), n}, srv{sr.data(), n},
                gxv{gx.data(), n}, grv{gr.data(), n}, fv{flux.data(), nf};
            for (int i = 0; i < n; ++i)
                wv.set(i, primitive(uv.get(i), settings.gas));
            for (int i = 0; i < n; ++i) {
                slopes(i, wv, sxv, srv, mesh.cells.data(), mesh.faces.data(), settings);
                transport_gradients(i, wv, gxv, grv, mesh.cells.data(), mesh.faces.data(),
                                    settings);
            }
            for (int f = 0; f < nf; ++f) {
                bool corrected = false;
                fv.set(f, fv_flux(f, wv, sxv, srv, gxv, grv, mesh.cells.data(), mesh.faces.data(),
                                  settings, corrected)
                              .value);
            }
            double inlet = 0, outlet = 0, wall = 0;
            for (int i = 0; i < nf; ++i) {
                auto f = mesh.faces[i];
                double mass = fv.get(i)[0] * f.area;
                if (f.boundary == Boundary::Inlet)
                    inlet += mass;
                if (f.boundary == Boundary::Outlet)
                    outlet += mass;
                if (f.boundary == Boundary::Wall)
                    wall += mass;
            }
            std::cout << suffix << " boundary mass flux: inlet=" << inlet << " outlet=" << outlet
                      << " wall=" << wall << " relative mismatch="
                      << std::abs(inlet - outlet) / std::max(std::abs(inlet), std::abs(outlet))
                      << '\n';
            std::ofstream out(prefix + suffix + ".csv");
            out << std::setprecision(17) << "i,j,x,r,rho,u,v,p,mass,axial,radial,energy\n";
            std::array<double, 4> sum{}, integral{};
            for (int i = 0; i < n; ++i) {
                auto rhs = fv_residual(i, fv, wv, gxv, grv, mesh.cells.data(), mesh.faces.data(),
                                       settings);
                auto q = wv.get(i);
                auto c = mesh.cells[i];
                out << i % mesh.nx << ',' << i / mesh.nx << ',' << c.x << ',' << c.r;
                for (int k = 0; k < 4; ++k)
                    out << ',' << q[k];
                for (int k = 0; k < 4; ++k) {
                    out << ',' << rhs[k];
                    sum[k] += rhs[k] * rhs[k] / n;
                    integral[k] += rhs[k] * c.volume;
                }
                out << '\n';
            }
            std::cout << suffix << " RMS:";
            for (double x : sum)
                std::cout << ' ' << std::sqrt(x);
            std::cout << " integrated RHS:";
            for (double x : integral)
                std::cout << ' ' << x;
            std::cout << '\n';
        };
        diagnose("-initial");
        for (int i = 0; i < steps; ++i) {
            solver->step();
            if (i == 0 || (i + 1) % 20 == 0) {
                auto s = solver->stats();
                log << s.iterations << ',' << s.time;
                for (double r : s.residual)
                    log << ',' << r;
                log << '\n';
            }
        }
        diagnose("-final");
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}

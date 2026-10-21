#include "astraflow/geometry/nozzle.hpp"
#include "astraflow/numerics/finite_volume.hpp"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <iostream>
#include <vector>
using namespace astraflow;
using Catch::Approx;
static void couette(bool gpu = false, bool fp64 = true) {
    auto m = rectangular_mesh(4, 24, 1, 1, true, false, true);
    Settings s;
    s.no_slip = true;
    s.wall_temperature = 1;
    s.upper_wall_speed = 1;
    s.gas.viscosity = 0.05;
    std::unique_ptr<Solver> solver;
#ifdef ASTRAFLOW_HAS_CUDA
    if (gpu)
        solver = make_cuda_solver(m, s, fp64);
#endif
    if (!solver)
        solver = make_cpu_solver(m, s, fp64);
    std::vector<State<double>> initial;
    for (auto c : m.cells)
        initial.push_back({{1, 0.8 * c.r, 0, 1}});
    solver->initialize(initial);
    while (solver->stats().time < 15)
        solver->step(15 - solver->stats().time);
    auto u = solver->state();
    int n = int(m.cells.size());
    double velocity_error = 0, temperature_error = 0;
    double conductivity = s.gas.viscosity * s.gas.gamma / ((s.gas.gamma - 1) * s.gas.prandtl);
    for (int i = 0; i < n; ++i) {
        State<double> q;
        for (int k = 0; k < 4; ++k)
            q[k] = u[k * n + i];
        auto w = primitive(q, s.gas);
        double y = m.cells[i].r, T = 1 + s.gas.viscosity / (2 * conductivity) * y * (1 - y);
        velocity_error += std::abs(w[1] - y) / n;
        temperature_error += std::abs(w[3] / w[0] - T) / n;
        REQUIRE(physical(w, s.gas));
    }
    std::cout << "Couette " << (gpu ? "CUDA" : "CPU") << " FP" << (fp64 ? 64 : 32)
              << " 4x24 t=15 steps=" << solver->stats().iterations
              << " velocity_L1=" << velocity_error << " temperature_L1=" << temperature_error
              << '\n';
    REQUIRE(velocity_error < 0.002);
    REQUIRE(temperature_error < 0.001);
}
TEST_CASE("Compressible Couette velocity and temperature", "[viscous]") { couette(); }
// With temperature-dependent viscosity the exact steady Couette solution still exists: the
// energy equation gives T(u) = Tw + Pr u (U - u) / (2 cp) for any mu(T) because k = mu cp / Pr,
// and the momentum equation mu(T(u)) du/dy = tau (constant) gives y(u) by quadrature.
static void variable_viscosity_couette(ViscosityModel model) {
    auto m = rectangular_mesh(4, 32, 1, 1, true, false, true);
    Settings s;
    s.no_slip = true;
    s.wall_temperature = 1;
    s.upper_wall_speed = 2;
    s.gas.viscosity = 0.05;
    s.gas.viscosity_model = model;
    s.gas.reference_temperature = 1;
    s.gas.sutherland_temperature = 0.5;
    s.gas.viscosity_exponent = 1;
    double cp = s.gas.gamma * s.gas.gas_constant / (s.gas.gamma - 1), U = s.gas.gas_constant;
    U = s.upper_wall_speed;
    auto exact_T = [&](double u) { return 1 + s.gas.prandtl * u * (U - u) / (2 * cp); };
    auto mu = [&](double u) { return viscosity(exact_T(u), s.gas); };
    // y(u) = (1/tau) int_0^u mu ds with tau fixed by y(U) = 1 (composite Simpson, 20000 panels).
    constexpr int panels = 20000;
    std::vector<double> integral(panels + 1, 0);
    double h = U / panels;
    for (int i = 1; i <= panels; ++i) {
        double a = (i - 1) * h, b = i * h, c = 0.5 * (a + b);
        integral[i] = integral[i - 1] + h / 6 * (mu(a) + 4 * mu(c) + mu(b));
    }
    double tau = integral[panels];
    auto exact_u = [&](double y) {
        int lo = 0, hi = panels;
        while (hi - lo > 1) {
            int mid = (lo + hi) / 2;
            (integral[mid] / tau <= y ? lo : hi) = mid;
        }
        double y0 = integral[lo] / tau, y1 = integral[hi] / tau;
        return (lo + (y - y0) / (y1 - y0)) * h;
    };
    REQUIRE(mu(U / 2) > 1.05 * mu(0)); // the model must actually vary across the channel
    auto solver = make_cpu_solver(m, s, true);
    std::vector<State<double>> initial;
    for (auto c : m.cells)
        initial.push_back({{1, 0.8 * U * c.r, 0, 1}});
    solver->initialize(initial);
    while (solver->stats().time < 20)
        solver->step(20 - solver->stats().time);
    auto u = solver->state();
    int n = int(m.cells.size());
    double velocity_error = 0, temperature_error = 0;
    for (int i = 0; i < n; ++i) {
        State<double> q;
        for (int k = 0; k < 4; ++k)
            q[k] = u[k * n + i];
        auto w = primitive(q, s.gas);
        double y = m.cells[i].r, ue = exact_u(y);
        velocity_error += std::abs(w[1] - ue) / (U * n);
        temperature_error += std::abs(w[3] / w[0] - exact_T(ue)) / n;
    }
    std::cout << "Variable-viscosity Couette "
              << (model == ViscosityModel::Sutherland ? "Sutherland" : "power-law")
              << " mu(centre)/mu(wall)=" << mu(U / 2) / mu(0) << " velocity_L1/U=" << velocity_error
              << " temperature_L1=" << temperature_error << '\n';
    REQUIRE(velocity_error < 0.002);
    REQUIRE(temperature_error < 0.001);
}
TEST_CASE("Couette flow with Sutherland viscosity matches the exact quadrature solution",
          "[viscous]") {
    variable_viscosity_couette(ViscosityModel::Sutherland);
}
TEST_CASE("Couette flow with power-law viscosity matches the exact quadrature solution",
          "[viscous]") {
    variable_viscosity_couette(ViscosityModel::PowerLaw);
}
TEST_CASE("Constant viscosity model reproduces the constant-viscosity arithmetic exactly",
          "[unit][viscous]") {
    Gas g;
    g.viscosity = 3e-5;
    for (double t : {0.1, 1.0, 7.5})
        REQUIRE(viscosity(t, g) == g.viscosity);
    g.viscosity_model = ViscosityModel::Sutherland;
    g.reference_temperature = 273.15;
    g.sutherland_temperature = 110.4;
    REQUIRE(viscosity(273.15, g) == Approx(3e-5));
    REQUIRE(viscosity(546.3, g) ==
            Approx(3e-5 * std::pow(2, 1.5) * (273.15 + 110.4) / (546.3 + 110.4)));
    g.viscosity_model = ViscosityModel::PowerLaw;
    g.viscosity_exponent = 0.7;
    REQUIRE(viscosity(546.3, g) == Approx(3e-5 * std::pow(2, 0.7)));
    Settings bad;
    bad.gas.viscosity_model = ViscosityModel::Sutherland;
    REQUIRE_THROWS(validate_settings(bad)); // reference temperature missing
}
#ifdef ASTRAFLOW_HAS_CUDA
TEST_CASE("CUDA long-time Couette analytical verification", "[cuda]") {
    couette(true, true);
    couette(true, false);
}
#endif
TEST_CASE("Viscous adiabatic rocket nozzle remains physical", "[viscous]") {
    Nozzle g;
    g.nx = 48;
    g.nr = 12;
    auto m = nozzle_mesh(g);
    Settings s;
    s.no_slip = true;
    s.gas.viscosity = 1e-4;
    auto solver = make_cpu_solver(m, s);
    solver->initialize(nozzle_initial_state(m, g, s));
    for (int i = 0; i < 1000; ++i)
        solver->step();
    auto u = solver->state();
    int n = int(m.cells.size());
    for (int i = 0; i < n; ++i) {
        State<double> q;
        for (int k = 0; k < 4; ++k)
            q[k] = u[k * n + i];
        REQUIRE(physical(primitive(q, s.gas), s.gas));
    }
    std::cout << "Viscous nozzle t=" << solver->stats().time
              << " reconstruction_fallbacks=" << solver->stats().reconstruction_fallbacks << '\n';
}
#ifdef ASTRAFLOW_HAS_CUDA
TEST_CASE("Viscous CUDA CPU parity in both precisions", "[cuda]") {
    Nozzle g;
    g.nx = 32;
    g.nr = 8;
    auto m = nozzle_mesh(g);
    Settings s;
    s.no_slip = true;
    s.gas.viscosity = 1e-4;
    for (bool fp64 : {true, false}) {
        auto cpu = make_cpu_solver(m, s, fp64);
        auto gpu = make_cuda_solver(m, s, fp64);
        auto w = nozzle_initial_state(m, g, s);
        cpu->initialize(w);
        gpu->initialize(w);
        for (int i = 0; i < 100; ++i) {
            cpu->step();
            gpu->step();
        }
        auto a = cpu->state(), b = gpu->state();
        double error = 0;
        for (std::size_t i = 0; i < a.size(); ++i)
            error = std::max(error, std::abs(a[i] - b[i]));
        std::cout << "Viscous FP" << (fp64 ? 64 : 32) << " max parity=" << error << '\n';
        REQUIRE(error < (fp64 ? 1e-10 : 5e-5));
    }
}
TEST_CASE("Sutherland viscosity CUDA CPU parity on an isothermal-wall nozzle", "[cuda]") {
    Nozzle g;
    g.nx = 32;
    g.nr = 8;
    auto m = nozzle_mesh(g);
    Settings s;
    s.no_slip = true;
    s.wall_temperature = 0.6;
    s.gas.viscosity = 1e-4;
    s.gas.viscosity_model = ViscosityModel::Sutherland;
    s.gas.reference_temperature = 1;
    s.gas.sutherland_temperature = 0.4;
    auto cpu = make_cpu_solver(m, s, true);
    auto gpu = make_cuda_solver(m, s, true);
    auto w = nozzle_initial_state(m, g, s);
    cpu->initialize(w);
    gpu->initialize(w);
    for (int i = 0; i < 100; ++i) {
        cpu->step();
        gpu->step();
    }
    auto a = cpu->state(), b = gpu->state();
    double error = 0;
    for (std::size_t i = 0; i < a.size(); ++i)
        error = std::max(error, std::abs(a[i] - b[i]));
    std::cout << "Sutherland FP64 max parity=" << error << '\n';
    REQUIRE(error < 1e-10);
}
TEST_CASE("CUDA profiling preserves solution and reports stage intervals", "[cuda]") {
    Nozzle g;
    g.nx = 16;
    g.nr = 8;
    auto mesh = nozzle_mesh(g);
    Settings s;
    s.no_slip = true;
    s.gas.viscosity = 1e-4;
    auto plain = make_cuda_solver(mesh, s, true);
    s.profile = true;
    auto profiled = make_cuda_solver(mesh, s, true);
    auto w = nozzle_initial_state(mesh, g, s);
    plain->initialize(w);
    profiled->initialize(w);
    for (int i = 0; i < 5; ++i) {
        plain->step();
        profiled->step();
    }
    REQUIRE(plain->state() == profiled->state());
    for (double ms : profiled->stats().stage_ms) {
        REQUIRE(std::isfinite(ms));
        REQUIRE(ms > 0);
    }
}
#endif

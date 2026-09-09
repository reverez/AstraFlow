#include "astraflow/geometry/nozzle.hpp"
#include "astraflow/numerics/finite_volume.hpp"
#include <catch2/catch_test_macros.hpp>
#include <iostream>
using namespace astraflow;
TEST_CASE("Compressible Couette velocity and temperature", "[viscous]") {
    auto m = rectangular_mesh(4, 24, 1, 1, true, false, true);
    Settings s;
    s.no_slip = true;
    s.wall_temperature = 1;
    s.upper_wall_speed = 1;
    s.gas.viscosity = 0.05;
    auto solver = make_cpu_solver(m, s);
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
    std::cout << "Couette 4x24 t=15 steps=" << solver->stats().iterations
              << " velocity_L1=" << velocity_error << " temperature_L1=" << temperature_error
              << '\n';
    REQUIRE(velocity_error < 0.002);
    REQUIRE(temperature_error < 0.001);
}
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
#endif

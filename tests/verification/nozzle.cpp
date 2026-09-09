#include "astraflow/geometry/nozzle.hpp"
#include "astraflow/numerics/finite_volume.hpp"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <iostream>
using namespace astraflow;
using Catch::Approx;
TEST_CASE("Axisymmetric metric closure and boundary invariants", "[unit]") {
    Nozzle g;
    g.nx = 32;
    g.nr = 8;
    auto m = nozzle_mesh(g);
    for (auto c : m.cells) {
        double x = 0, r = 0;
        for (int d = 0; d < 4; ++d) {
            auto f = m.faces[c.faces[d]];
            double sign = d % 2 ? 1 : -1;
            x += sign * f.area * f.nx;
            r += sign * f.area * f.nr;
        }
        REQUIRE(std::abs(x) < 1e-14);
        REQUIRE(r == Approx(c.planar_area).margin(1e-14));
        REQUIRE(c.r > 0);
    }
    REQUIRE(g.derivative(g.chamber_length) == Approx(0));
    REQUIRE(g.derivative(g.chamber_length + g.contraction_length) == Approx(0));
    g.throat_radius = -1;
    REQUIRE_THROWS(nozzle_mesh(g));
    Settings s;
    Face f;
    f.boundary = Boundary::Inlet;
    auto inlet = boundary_state(State<double>{{1, 0.2, 0, 1}}, f, s);
    double t = inlet[3] / inlet[0];
    REQUIRE(t + (s.gas.gamma - 1) / (2 * s.gas.gamma) * inlet[1] * inlet[1] == Approx(s.t0));
    f.boundary = Boundary::Outlet;
    auto sup = State<double>{{1, 3, 0, 1}};
    auto out = boundary_state(sup, f, s);
    for (int k = 0; k < 4; ++k)
        REQUIRE(out[k] == sup[k]);
}
TEST_CASE("Axis remains finite and uniform pressure balances geometry", "[verification]") {
    Nozzle g;
    g.nx = 32;
    g.nr = 8;
    auto m = nozzle_mesh(g);
    Settings s;
    s.back_pressure = 1;
    auto solver = make_cpu_solver(m, s);
    solver->initialize(std::vector<State<double>>(m.cells.size(), {{1, 0, 0, 1}}));
    auto before = solver->state();
    for (int i = 0; i < 50; ++i)
        solver->step();
    auto after = solver->state();
    double maximum = 0;
    for (std::size_t i = 0; i < before.size(); ++i)
        maximum = std::max(maximum, std::abs(after[i] - before[i]));
    std::cout << "Axis rest max drift=" << maximum << '\n';
    REQUIRE(maximum < 1e-11);
}
TEST_CASE("Inviscid nozzle choking and quasi-1D comparison", "[nozzle]") {
    Nozzle g;
    auto m = nozzle_mesh(g);
    Settings s;
    s.back_pressure = 0.07;
    auto solver = make_cpu_solver(m, s);
    solver->initialize(nozzle_initial_state(m, g, s));
    for (int i = 0; i < 4000; ++i)
        solver->step();
    auto u = solver->state();
    int n = int(m.cells.size());
    double minmass = 1e100, maxmass = 0, error = 0, throat = 0, near = 1e100;
    for (int i = 0; i < g.nx; ++i) {
        double mass = 0, mach = 0, area = 0;
        for (int j = 0; j < g.nr; ++j) {
            int k = j * g.nx + i;
            State<double> q;
            for (int d = 0; d < 4; ++d)
                q[d] = u[d * n + k];
            auto w = primitive(q, s.gas);
            REQUIRE(physical(w, s.gas));
            double a = m.cells[k].volume / m.cells[k].dx;
            mass += w[0] * w[1] * a;
            mach += w[1] / sound_speed(w, s.gas) * a;
            area += a;
        }
        mach /= area;
        minmass = std::min(minmass, mass);
        maxmass = std::max(maxmass, mass);
        double x = m.cells[i].x, R = g.radius(x),
               ref = isentropic_mach(R * R / (g.throat_radius * g.throat_radius), s.gas.gamma,
                                     x > g.chamber_length + g.contraction_length);
        error += std::abs(mach - ref) / g.nx;
        double distance = std::abs(x - g.chamber_length - g.contraction_length);
        if (distance < near) {
            near = distance;
            throat = mach;
        }
    }
    double spread = (maxmass - minmass) / maxmass;
    std::cout << "Nozzle 96x12 steps=4000 t=" << solver->stats().time << " throat_M=" << throat
              << " Mach_L1=" << error << " mass_spread=" << spread << '\n';
    REQUIRE(throat > 0.85);
    REQUIRE(throat < 1.15);
    REQUIRE(error < 0.12);
    REQUIRE(spread < 0.05);
}
#ifdef ASTRAFLOW_HAS_CUDA
TEST_CASE("CUDA axisymmetric nozzle parity", "[cuda]") {
    Nozzle g;
    g.nx = 32;
    g.nr = 8;
    auto m = nozzle_mesh(g);
    Settings s;
    auto cpu = make_cpu_solver(m, s);
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
    std::cout << "Axisymmetric FP64 max parity=" << error << '\n';
    REQUIRE(error < 1e-10);
}
#endif

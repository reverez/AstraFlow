#include "astraflow/core/solver.hpp"
#include <catch2/catch_test_macros.hpp>
#include <iostream>
#include <numbers>
using namespace astraflow;
TEST_CASE("Smooth 2D entropy wave converges under grid refinement", "[refinement]") {
    double errors[3]{};
    int level = 0;
    for (int nx : {20, 40, 80}) {
        auto mesh = rectangular_mesh(nx, nx / 2, 1, 1, true, true);
        Settings settings;
        auto solver = make_cpu_solver(mesh, settings);
        std::vector<State<double>> initial;
        for (auto c : mesh.cells)
            initial.push_back(
                {{1 + 0.2 * std::sin(2 * std::numbers::pi * (c.x + c.r)), 0.4, 0.2, 1}});
        solver->initialize(initial);
        while (solver->stats().time < 0.25)
            solver->step(0.25 - solver->stats().time);
        auto u = solver->state();
        double l1 = 0, l2 = 0, maximum = 0, mass = 0;
        int n = int(mesh.cells.size());
        for (int i = 0; i < n; ++i) {
            auto c = mesh.cells[i];
            double exact = 1 + 0.2 * std::sin(2 * std::numbers::pi * (c.x + c.r - 0.6 * 0.25)),
                   difference = std::abs(u[i] - exact);
            l1 += difference / n;
            l2 += difference * difference / n;
            maximum = std::max(maximum, difference);
            mass += u[i] * c.volume;
        }
        errors[level++] = l1;
        std::cout << "Refinement entropy " << nx << 'x' << nx / 2 << " L1=" << l1
                  << " L2=" << std::sqrt(l2) << " Linf=" << maximum
                  << " mass_error=" << std::abs(mass - 1) << '\n';
        REQUIRE(std::abs(mass - 1) < 1e-12);
    }
    double p1 = std::log2(errors[0] / errors[1]), p2 = std::log2(errors[1] / errors[2]);
    std::cout << "Entropy observed orders=" << p1 << ", " << p2 << '\n';
    REQUIRE(p1 > 1.5);
    REQUIRE(p2 > 1.5);
    REQUIRE(errors[2] < 0.002);
}
TEST_CASE("Compressible Couette temperature converges under refinement", "[refinement]") {
    double errors[3]{};
    int level = 0;
    for (int nr : {12, 24, 48}) {
        auto mesh = rectangular_mesh(4, nr, 1, 1, true, false, true);
        Settings s;
        s.no_slip = true;
        s.wall_temperature = 1;
        s.upper_wall_speed = 1;
        s.gas.viscosity = 0.05;
        auto solver = make_cpu_solver(mesh, s);
        std::vector<State<double>> initial;
        for (auto c : mesh.cells)
            initial.push_back({{1, 0.8 * c.r, 0, 1}});
        solver->initialize(initial);
        while (solver->stats().time < 20)
            solver->step(20 - solver->stats().time);
        auto u = solver->state();
        int n = int(mesh.cells.size());
        double l1 = 0, l2 = 0, maximum = 0, velocity = 0,
               k = s.gas.viscosity * s.gas.gamma / ((s.gas.gamma - 1) * s.gas.prandtl);
        for (int i = 0; i < n; ++i) {
            State<double> q;
            for (int d = 0; d < 4; ++d)
                q[d] = u[d * n + i];
            auto w = primitive(q, s.gas);
            double y = mesh.cells[i].r, exact = 1 + s.gas.viscosity / (2 * k) * y * (1 - y),
                   d = std::abs(w[3] / w[0] - exact);
            l1 += d / n;
            l2 += d * d / n;
            maximum = std::max(maximum, d);
            velocity += std::abs(w[1] - y) / n;
        }
        errors[level++] = l1;
        std::cout << "Refinement Couette 4x" << nr << " temperature_L1=" << l1
                  << " L2=" << std::sqrt(l2) << " Linf=" << maximum << " velocity_L1=" << velocity
                  << '\n';
        REQUIRE(velocity < 1e-4);
    }
    double p1 = std::log2(errors[0] / errors[1]), p2 = std::log2(errors[1] / errors[2]);
    std::cout << "Couette observed orders=" << p1 << ", " << p2 << '\n';
    REQUIRE(p1 > 1.7);
    REQUIRE(p2 > 1.7);
    REQUIRE(errors[2] < 3e-5);
}

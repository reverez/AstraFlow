#include "astraflow/core/solver.hpp"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <iostream>
using namespace astraflow;
using Catch::Approx;
TEST_CASE("Rectangular mesh geometry and face closure", "[unit]") {
    auto m = rectangular_mesh(8, 4, 2, 1);
    double volume = 0;
    for (auto c : m.cells) {
        volume += c.volume;
        double x = 0, r = 0;
        for (int d = 0; d < 4; ++d) {
            auto f = m.faces[c.faces[d]];
            double sign = d % 2 ? 1 : -1;
            x += sign * f.area * f.nx;
            r += sign * f.area * f.nr;
        }
        REQUIRE(std::abs(x) < 1e-14);
        REQUIRE(std::abs(r) < 1e-14);
    }
    REQUIRE(volume == Approx(2));
    REQUIRE_THROWS(rectangular_mesh(0, 2));
    REQUIRE_THROWS(rectangular_mesh(8, 8, -1));
    m.cells[0].volume = -1;
    REQUIRE_THROWS(m.validate());
}
static void uniform(bool gpu) {
    auto m = rectangular_mesh(16, 8, 1, 1, true, true);
    Settings s;
    std::unique_ptr<Solver> solver;
#ifdef ASTRAFLOW_HAS_CUDA
    if (gpu)
        solver = make_cuda_solver(m, s, true);
#endif
    if (!solver)
        solver = make_cpu_solver(m, s);
    solver->initialize(std::vector<State<double>>(m.cells.size(), {{1, 0.4, -0.2, 1}}));
    auto before = solver->state();
    for (int i = 0; i < 20; ++i)
        solver->step();
    auto after = solver->state();
    for (std::size_t i = 0; i < after.size(); ++i)
        REQUIRE(std::abs(after[i] - before[i]) < 1e-13);
}
TEST_CASE("CPU 2D uniform flow", "[verification]") { uniform(false); }
#ifdef ASTRAFLOW_HAS_CUDA
TEST_CASE("CUDA 2D uniform flow", "[cuda]") { uniform(true); }
TEST_CASE("2D conservative evolution and CPU GPU parity", "[cuda]") {
    auto m = rectangular_mesh(32, 16, 1, 1, true, true);
    Settings s;
    auto cpu = make_cpu_solver(m, s);
    auto gpu = make_cuda_solver(m, s, true);
    std::vector<State<double>> w;
    for (auto c : m.cells)
        w.push_back(
            {{1 + 0.1 * std::sin(6.283185307179586 * c.x) * std::cos(6.283185307179586 * c.r), 0.4,
              0.2, 1}});
    cpu->initialize(w);
    gpu->initialize(w);
    auto initial = cpu->state();
    for (int i = 0; i < 30; ++i) {
        cpu->step();
        gpu->step();
    }
    auto a = cpu->state(), b = gpu->state();
    int n = int(m.cells.size());
    double maximum = 0;
    for (int k = 0; k < 4; ++k) {
        double before = 0, after = 0;
        for (int i = 0; i < n; ++i) {
            before += initial[k * n + i] * m.cells[i].volume;
            after += a[k * n + i] * m.cells[i].volume;
            maximum = std::max(maximum, std::abs(a[k * n + i] - b[k * n + i]));
        }
        REQUIRE(std::abs(after - before) < 1e-12);
    }
    std::cout << "2D FP64 max parity=" << maximum << '\n';
    REQUIRE(maximum < 2e-11);
}
#endif

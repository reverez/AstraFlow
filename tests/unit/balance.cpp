#include "astraflow/analysis/balance.hpp"
#include "astraflow/core/simulation.hpp"
#include "astraflow/numerics/flux.hpp"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
using namespace astraflow;
TEST_CASE("FV balance includes pressure wall and cylindrical source", "[unit][balance]") {
    Nozzle shape;
    shape.nx = 16;
    shape.nr = 8;
    auto mesh = nozzle_mesh(shape);
    Settings s;
    s.p0 = s.back_pressure = 1;
    s.t0 = 1;
    int n = int(mesh.cells.size());
    std::vector<double> u(4 * n);
    auto q = conservative(State<double>{{1, 0, 0, 1}}, s.gas);
    for (int i = 0; i < n; ++i)
        for (int k = 0; k < 4; ++k)
            u[k * n + i] = q[k];
    auto b = finite_volume_balance(mesh, u, s);
    REQUIRE(b.global["local_decomposition_max_error"].get<double>() < 1e-12);
    for (int k = 0; k < 4; ++k) {
        REQUIRE(b.global["rms"][k].get<double>() < 1e-11);
        REQUIRE(std::abs(b.global["assembly_identity_error"][k].get<double>()) < 1e-12);
    }
    REQUIRE(b.global["pressure_source"][2].get<double>() > 0);
    REQUIRE(b.global["boundary_fluxes"]["wall"]["outward_convective"][2].get<double>() > 0);
}
TEST_CASE("Viscous energy and all four balances decompose consistently", "[unit][balance]") {
    auto m = rectangular_mesh(12, 8, 1, 1, false, false, true);
    Settings s;
    s.no_slip = true;
    s.upper_wall_speed = 1;
    s.wall_temperature = 1;
    s.gas.viscosity = .05;
    int n = int(m.cells.size());
    std::vector<double> u(4 * n);
    for (int i = 0; i < n; ++i) {
        auto c = m.cells[i];
        auto q = conservative(
            State<double>{{1 + .01 * c.x, .8 * c.r, .01 * std::sin(c.x), 1 + .03 * c.r}}, s.gas);
        for (int k = 0; k < 4; ++k)
            u[k * n + i] = q[k];
    }
    auto b = finite_volume_balance(m, u, s);
    for (int k = 0; k < 4; ++k)
        REQUIRE(std::abs(b.global["assembly_identity_error"][k].get<double>()) < 1e-12);
    REQUIRE(b.global["local_decomposition_max_error"].get<double>() < 1e-12);
    REQUIRE(
        std::abs(b.global["boundary_fluxes"]["wall"]["outward_viscous_thermal"][3].get<double>()) >
        1e-4);
}
TEST_CASE("Conservative prolongation preserves each parent's four integrals", "[unit][balance]") {
    Nozzle c;
    c.nx = 16;
    c.nr = 8;
    auto f = c;
    f.nx *= 2;
    f.nr *= 2;
    auto coarse = nozzle_mesh(c), fine = nozzle_mesh(f);
    Gas gas;
    int n = int(coarse.cells.size());
    std::vector<double> u(4 * n);
    for (int i = 0; i < n; ++i) {
        auto cell = coarse.cells[i];
        auto q = conservative(
            State<double>{{1 + .05 * cell.x, 1 + .1 * cell.r, .02, 1 + .03 * cell.x}}, gas);
        for (int k = 0; k < 4; ++k)
            u[k * n + i] = q[k];
    }
    nlohmann::json audit;
    auto out = conservative_prolong(coarse, fine, u, gas, audit);
    REQUIRE(audit["maximum_parent_relative_integral_error"].get<double>() < 1e-14);
    REQUIRE(audit["maximum_volume_ratio_change"].get<double>() > 0);
    for (auto q : out)
        REQUIRE(physical(primitive(q, gas), gas));
    for (int j = 0; j < c.nr; ++j)
        for (int i = 0; i < c.nx; ++i)
            for (int k = 0; k < 4; ++k) {
                double sum = 0;
                for (int b = 0; b < 2; ++b)
                    for (int a = 0; a < 2; ++a) {
                        int id = (2 * j + b) * f.nx + 2 * i + a;
                        sum += out[id][k] * fine.cells[id].volume;
                    }
                REQUIRE(sum ==
                        Catch::Approx(u[k * n + j * c.nx + i] * coarse.cells[j * c.nx + i].volume)
                            .epsilon(1e-13));
            }
    REQUIRE(limit(1., 2., Limiter::FirstOrderDiagnostic) == 0);
    REQUIRE(limit(1., 2., Limiter::MC) == 1.5);
    REQUIRE_THROWS(Config::parse({{"numerics", {{"limiter", "first_order"}}}}));
}
TEST_CASE("Snapshot balance is the production time derivative", "[unit][balance]") {
    Nozzle geometry;
    geometry.nx = 16;
    geometry.nr = 8;
    auto mesh = nozzle_mesh(geometry);
    Settings s;
    s.no_slip = true;
    s.gas.viscosity = .0001;
    auto initial = nozzle_initial_state(mesh, geometry, s);
    for (auto limiter : {Limiter::MC, Limiter::VanLeer, Limiter::FirstOrderDiagnostic}) {
        s.limiter = limiter;
        auto solver = make_cpu_solver(mesh, s, true);
        solver->initialize(initial);
        auto before = solver->state();
        auto balance = finite_volume_balance(mesh, before, s);
        double dt = 1e-9;
        solver->step(dt);
        auto after = solver->state();
        for (std::size_t i = 0; i < mesh.cells.size(); ++i)
            for (int k = 0; k < 4; ++k) {
                auto index = k * mesh.cells.size() + i;
                double derivative = (after[index] - before[index]) / dt;
                REQUIRE(std::abs(derivative - balance.total[i][k]) <
                        2e-5 * std::max(1.0, std::abs(balance.total[i][k])));
            }
    }
}
TEST_CASE("Explicit physical initialization preserves all conservative components",
          "[unit][balance]") {
    Config c;
    c.backend = "cpu";
    c.precision = "double";
    c.geometry.nx = 8;
    c.geometry.nr = 4;
    auto m = nozzle_mesh(c.geometry);
    auto w = nozzle_initial_state(m, c.geometry, c.settings);
    for (auto &q : w) {
        q[1] *= .9;
        q[2] += .01;
    }
    Simulation sim(c, w);
    auto u = sim.state();
    for (std::size_t i = 0; i < w.size(); ++i) {
        auto expected = conservative(w[i], c.settings.gas);
        for (int k = 0; k < 4; ++k)
            REQUIRE(u[k * w.size() + i] == Catch::Approx(expected[k]).epsilon(1e-14));
    }
    REQUIRE(sim.stats().iterations == 0);
    REQUIRE_FALSE(sim.convergence_report()["converged"].get<bool>());
    w.pop_back();
    REQUIRE_THROWS(Simulation(c, w));
}

#include "astraflow/analysis/balance.hpp"
#include "astraflow/geometry/nozzle.hpp"
#include <algorithm>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <iostream>
using namespace astraflow;
using Catch::Approx;
namespace {
struct Case {
    Mesh mesh;
    Settings settings;
    std::vector<State<double>> initial;
};
Case viscous_nozzle(int nx = 24, int nr = 8) {
    Nozzle g;
    g.nx = nx;
    g.nr = nr;
    Case c{nozzle_mesh(g), {}, {}};
    c.settings.no_slip = true;
    c.settings.gas.viscosity = 1e-4;
    c.initial = nozzle_initial_state(c.mesh, g, c.settings);
    return c;
}
std::array<double, 2> host_boundary_mass(const Mesh &m, const std::vector<double> &u, Settings s) {
    auto g = finite_volume_balance(m, u, s).global;
    const auto &flux = g.at("boundary_fluxes");
    return {-flux.at("inlet").at("outward_total")[0].get<double>(),
            flux.at("outlet").at("outward_total")[0].get<double>()};
}
} // namespace
TEST_CASE("CPU reference fields and statistics are independent of thread count",
          "[unit][parallel]") {
    auto c = viscous_nozzle(64, 40); // above the serial threshold, so loops really fork
    int workers = std::min(8, cpu_hardware_threads());
    std::cout << "OpenMP workers used: " << workers << " of " << cpu_hardware_threads() << '\n';
    set_cpu_threads(1);
    auto single = make_cpu_solver(c.mesh, c.settings, true);
    single->initialize(c.initial);
    set_cpu_threads(workers);
    auto parallel = make_cpu_solver(c.mesh, c.settings, true);
    parallel->initialize(c.initial);
    for (int i = 0; i < 30; ++i) {
        single->step();
        parallel->step();
    }
    REQUIRE(single->state() == parallel->state());
    REQUIRE(single->stats().residual == parallel->stats().residual);
    REQUIRE(single->stats().boundary_mass_flux == parallel->stats().boundary_mass_flux);
    REQUIRE(single->stats().flux_fallbacks == parallel->stats().flux_fallbacks);
    REQUIRE(parallel->stats().threads == workers);
    REQUIRE_THROWS(set_cpu_threads(0));
    set_cpu_threads(1);
}
TEST_CASE("Pseudo-mode boundary mass fluxes equal the host balance at the completed iterate",
          "[unit][parallel]") {
    auto c = viscous_nozzle();
    SteadySettings steady;
    steady.enabled = true;
    auto solver = make_cpu_solver(c.mesh, c.settings, true, steady);
    solver->initialize(c.initial);
    for (int i = 0; i < 25; ++i)
        solver->step();
    auto exact = host_boundary_mass(c.mesh, solver->state(), c.settings);
    auto measured = solver->stats().boundary_mass_flux;
    for (int k = 0; k < 2; ++k)
        REQUIRE(measured[k] == Approx(exact[k]).margin(1e-14));
    REQUIRE(boundary_mass_mismatch(measured) ==
            Approx(std::abs(exact[1] - exact[0]) / std::max(std::abs(exact[0]), std::abs(exact[1])))
                .margin(1e-14));
    REQUIRE(boundary_mass_mismatch({0, 0}) == 0);
    REQUIRE(boundary_mass_mismatch({1, 1.5}) == Approx(1.0 / 3));
    REQUIRE(boundary_mass_mismatch({-2, -1}) == Approx(0.5));
}
TEST_CASE("Physical-mode boundary mass fluxes come from the second stage, not the new state",
          "[unit][parallel]") {
    auto c = viscous_nozzle();
    auto solver = make_cpu_solver(c.mesh, c.settings, true);
    solver->initialize(c.initial);
    solver->step();
    auto stage = solver->stats().boundary_mass_flux;
    auto final_state = host_boundary_mass(c.mesh, solver->state(), c.settings);
    // The sampled values are finite, physically oriented (positive inflow/outflow for this
    // nozzle) and close to, but not necessarily identical with, the post-step evaluation.
    for (int k = 0; k < 2; ++k) {
        REQUIRE(std::isfinite(stage[k]));
        REQUIRE(stage[k] > 0);
        REQUIRE(stage[k] == Approx(final_state[k]).epsilon(1e-2));
    }
}
#ifdef ASTRAFLOW_HAS_CUDA
TEST_CASE("CUDA boundary mass fluxes match the CPU reference in both modes", "[cuda][parallel]") {
    auto c = viscous_nozzle(32, 8);
    for (bool pseudo : {false, true}) {
        SteadySettings steady;
        steady.enabled = pseudo;
        auto cpu = make_cpu_solver(c.mesh, c.settings, true, steady),
             gpu = make_cuda_solver(c.mesh, c.settings, true, steady);
        cpu->initialize(c.initial);
        gpu->initialize(c.initial);
        for (int i = 0; i < 40; ++i) {
            cpu->step();
            gpu->step();
        }
        auto a = cpu->stats().boundary_mass_flux, b = gpu->stats().boundary_mass_flux;
        std::cout << (pseudo ? "pseudo" : "physical") << " boundary mass flux CPU=" << a[0] << ','
                  << a[1] << " GPU=" << b[0] << ',' << b[1] << '\n';
        for (int k = 0; k < 2; ++k)
            REQUIRE(a[k] == Approx(b[k]).margin(1e-12));
        if (pseudo) {
            auto exact = host_boundary_mass(c.mesh, gpu->state(), c.settings);
            for (int k = 0; k < 2; ++k)
                REQUIRE(b[k] == Approx(exact[k]).margin(1e-12));
        }
        REQUIRE(gpu->stats().dt == Approx(cpu->stats().dt).margin(1e-15));
    }
}
#endif

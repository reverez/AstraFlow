#include "astraflow/analysis/wall.hpp"
#include "astraflow/analysis/balance.hpp"
#include "astraflow/geometry/nozzle.hpp"
#include "astraflow/io/output.hpp"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <iostream>
using namespace astraflow;
using Catch::Approx;
namespace {
// Steady compressible Couette flow with equal isothermal walls: exact shear mu U / H and exact
// wall heat flux mu U^2 / (2 H) into each wall (all dissipation leaves through the walls).
struct Couette {
    Mesh mesh = rectangular_mesh(4, 48, 1, 1, true, false, true);
    Settings s;
    std::vector<double> state;
    Couette() {
        s.no_slip = true;
        s.wall_temperature = 1;
        s.upper_wall_speed = 1;
        s.gas.viscosity = 0.05;
        auto solver = make_cpu_solver(mesh, s, true);
        std::vector<State<double>> initial;
        for (auto c : mesh.cells)
            initial.push_back({{1, 0.8 * c.r, 0, 1}});
        solver->initialize(initial);
        while (solver->stats().time < 20)
            solver->step(20 - solver->stats().time);
        state = solver->state();
    }
};
} // namespace
TEST_CASE("Couette wall shear and heat flux match the exact steady solution", "[unit][wall]") {
    Couette c;
    auto samples = wall_distribution(c.mesh, c.state, c.s);
    REQUIRE(samples.size() == 8); // four lower and four upper wall faces
    double shear = c.s.gas.viscosity * c.s.upper_wall_speed, heat = 0.5 * shear;
    int lower = 0, upper = 0;
    for (const auto &w : samples) {
        (w.upper ? upper : lower)++;
        // The stationary lower wall drags the fluid backwards; the moving upper wall pulls it
        // forward. At the discrete steady state energy conservation forces the wall conductive
        // flux to balance exactly half the discrete dissipation, so the heat flux reproduces
        // mu U^2 / (2H) to the steadiness of the solution, not to the wall-gradient truncation.
        REQUIRE(w.traction_x == Approx(w.upper ? shear : -shear).epsilon(2e-3));
        REQUIRE(w.traction_r == Approx(0).margin(1e-6 * shear)); // decaying transient v
        REQUIRE(w.heat_flux_into_wall == Approx(heat).epsilon(1e-3));
        REQUIRE(w.temperature == Approx(1).epsilon(1e-2));
        REQUIRE(w.length == Approx(0.25));
        REQUIRE(w.area == Approx(0.25));
    }
    REQUIRE(lower == 4);
    REQUIRE(upper == 4);
    auto summary = wall_summary(c.mesh, samples);
    std::cout << "Couette wall force=" << summary["wall_axial_viscous_force"]
              << " heat rate=" << summary["wall_heat_transfer_rate"]
              << " dissipation=" << shear * c.s.upper_wall_speed << '\n';
    // Both walls see the same traction with opposite outward normals, so net axial force is 0
    // and the total wall heat transfer equals the viscous dissipation mu U^2 / H.
    REQUIRE(summary["wall_axial_viscous_force"].get<double>() == Approx(0).margin(2e-3 * shear));
    REQUIRE(summary["wall_heat_transfer_rate"].get<double>() ==
            Approx(shear * c.s.upper_wall_speed).epsilon(2e-2));
    REQUIRE(summary["wall_peak_shear_stress"].get<double>() == Approx(shear).epsilon(2e-3));
}
TEST_CASE("Adiabatic nozzle wall has zero heat flux and a decelerating shear", "[unit][wall]") {
    Nozzle g;
    g.nx = 48;
    g.nr = 12;
    auto m = nozzle_mesh(g);
    Settings s;
    s.no_slip = true;
    s.gas.viscosity = 1e-4;
    auto solver = make_cpu_solver(m, s, true);
    solver->initialize(nozzle_initial_state(m, g, s));
    for (int i = 0; i < 500; ++i)
        solver->step();
    auto samples = wall_distribution(m, solver->state(), s);
    REQUIRE(samples.size() == std::size_t(g.nx));
    double force = 0;
    for (const auto &w : samples) {
        REQUIRE(w.upper);
        REQUIRE(w.heat_flux_into_wall == Approx(0).margin(1e-12));
        REQUIRE(std::isfinite(w.traction_x));
        force += w.traction_x * w.area;
    }
    // The fluid moves in +x, so the wall drags it backwards: traction on the fluid is negative.
    REQUIRE(force < 0);
    auto summary = wall_summary(m, samples);
    REQUIRE(summary["wall_heat_transfer_rate"].get<double>() == Approx(0).margin(1e-9));
    REQUIRE(summary["wall_axial_viscous_force"].get<double>() > 0);
    REQUIRE(summary["wall_axial_viscous_force"].get<double>() ==
            Approx(-2 * std::numbers::pi * force).epsilon(1e-12));
}
TEST_CASE("Isothermal nozzle wall heat transfer equals the wall energy flux of the balance",
          "[unit][wall]") {
    Nozzle g;
    g.nx = 32;
    g.nr = 8;
    auto m = nozzle_mesh(g);
    Settings s;
    s.no_slip = true;
    s.wall_temperature = 0.5;
    s.gas.viscosity = 1e-3;
    auto solver = make_cpu_solver(m, s, true);
    solver->initialize(nozzle_initial_state(m, g, s));
    for (int i = 0; i < 200; ++i)
        solver->step();
    auto state = solver->state();
    auto samples = wall_distribution(m, state, s);
    auto summary = wall_summary(m, samples);
    auto balance = finite_volume_balance(m, state, s).global;
    // Balance fluxes are outward and per radian; wall viscous/thermal energy flux leaving the
    // fluid is the same physical quantity as the integrated wall heat transfer rate.
    double balance_heat =
        -2 * std::numbers::pi *
        balance["boundary_fluxes"]["wall"]["outward_viscous_thermal"][3].get<double>();
    std::cout << "Isothermal nozzle wall heat rate=" << summary["wall_heat_transfer_rate"]
              << " balance=" << balance_heat << '\n';
    REQUIRE(summary["wall_heat_transfer_rate"].get<double>() ==
            Approx(balance_heat).epsilon(1e-10));
    REQUIRE(summary["wall_heat_transfer_rate"].get<double>() > 0); // hot gas cools to the wall
}
TEST_CASE("Isothermal walls are exactly impermeable", "[unit][wall]") {
    Nozzle g;
    g.nx = 32;
    g.nr = 8;
    auto m = nozzle_mesh(g);
    Settings s;
    s.no_slip = true;
    s.wall_temperature = 0.5;
    s.gas.viscosity = 1e-3;
    auto solver = make_cpu_solver(m, s, true);
    solver->initialize(nozzle_initial_state(m, g, s));
    for (int i = 0; i < 200; ++i)
        solver->step();
    auto balance = finite_volume_balance(m, solver->state(), s).global;
    // A cold wall changes the ghost density; the convective wall flux must still carry no mass.
    REQUIRE(balance["boundary_fluxes"]["wall"]["outward_convective"][0].get<double>() ==
            Approx(0).margin(1e-15));
    REQUIRE(balance["boundary_fluxes"]["wall"]["outward_convective"][3].get<double>() ==
            Approx(0).margin(1e-15));
    REQUIRE(balance["boundary_fluxes"]["wall"]["outward_viscous_thermal"][3].get<double>() != 0);
}
TEST_CASE("State export and import roundtrip exactly", "[unit][io]") {
    Config c;
    c.backend = "cpu";
    c.precision = "double";
    c.geometry.nx = 12;
    c.geometry.nr = 6;
    c.settings.no_slip = true;
    c.settings.gas.viscosity = 1e-4;
    Simulation sim(c);
    for (int i = 0; i < 5; ++i)
        sim.step();
    auto exported = export_state(sim);
    REQUIRE(exported["schema"] == 1);
    REQUIRE(exported["iterations"] == 5);
    auto initial = import_state(exported, c);
    Simulation restarted(c, initial);
    auto a = sim.state(), b = restarted.state();
    for (std::size_t i = 0; i < a.size(); ++i)
        REQUIRE(a[i] == Approx(b[i]).margin(1e-13));
    auto wrong = c;
    wrong.geometry.nx = 14;
    REQUIRE_THROWS(import_state(exported, wrong));
    auto tampered = exported;
    tampered["physical_conservative"][0][0] = -1;
    REQUIRE_THROWS(import_state(tampered, c));
    tampered = exported;
    tampered["schema"] = 2;
    REQUIRE_THROWS(import_state(tampered, c));
}

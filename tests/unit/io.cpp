#include "astraflow/io/output.hpp"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <numbers>
using namespace astraflow;
using Catch::Approx;
TEST_CASE("Configuration roundtrip and invalid fields", "[unit]") {
    Config c;
    c.backend = "cpu";
    c.precision = "double";
    auto round = Config::parse(c.json());
    REQUIRE(round.json() == c.json());
    auto j = c.json();
    j["gas"]["typo"] = 1;
    REQUIRE_THROWS(Config::parse(j));
    j = c.json();
    j["numerics"]["limiter"] = "none";
    REQUIRE_THROWS(Config::parse(j));
    j = c.json();
    j["gas"]["rho_floor"] = -1;
    REQUIRE_THROWS(Config::parse(j));
    j = c.json();
    j["runtime"]["max_iterations"] = 0;
    REQUIRE_THROWS(Config::parse(j));
}
TEST_CASE("Engineering exit integrals use actual annular area", "[unit]") {
    Nozzle g;
    g.nx = 16;
    g.nr = 8;
    auto m = nozzle_mesh(g);
    int n = int(m.cells.size());
    std::vector<double> u(4 * n);
    auto q = conservative(State<double>{{1, 2, 0, 3}}, Gas{});
    for (int k = 0; k < 4; ++k)
        for (int i = 0; i < n; ++i)
            u[k * n + i] = q[k];
    auto result = engineering(m, u, Gas{}, 1);
    double area = std::numbers::pi * g.exit_radius * g.exit_radius;
    REQUIRE(result["exit_area"].get<double>() == Approx(area));
    REQUIRE(result["mass_flow"].get<double>() == Approx(2 * area));
    REQUIRE(result["estimated_thrust"].get<double>() == Approx(6 * area));
    REQUIRE(result["specific_impulse"].get<double>() == Approx(3 / 9.80665));
    REQUIRE(result["exit_pressure"].get<double>() == Approx(3));
}
TEST_CASE("Dimensional reference scaling roundtrip", "[unit]") {
    Config c;
    c.backend = "cpu";
    c.precision = "double";
    c.geometry.nx = 16;
    c.geometry.nr = 8;
    c.settings.p0 = 2e6;
    c.settings.t0 = 2800;
    c.settings.gas.gas_constant = 287;
    c.settings.back_pressure = 101325;
    Simulation sim(c);
    auto expected = nozzle_initial_state(sim.mesh(), c.geometry, c.settings);
    auto u = sim.state();
    int n = int(expected.size());
    for (int i = 0; i < n; ++i) {
        State<double> q;
        for (int k = 0; k < 4; ++k)
            q[k] = u[k * n + i];
        auto w = primitive(q, c.settings.gas);
        for (int k = 0; k < 4; ++k)
            REQUIRE(w[k] == Approx(expected[i][k]).margin(1e-10));
    }
    sim.step();
    REQUIRE(sim.stats().time > 0);
    REQUIRE(sim.scales().pressure == 2e6);
}

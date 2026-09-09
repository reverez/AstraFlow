#include "../verification/exact_riemann.hpp"
#include "astraflow/core/euler1d.hpp"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <iostream>
#include <numeric>
using namespace astraflow;
using Catch::Approx;
TEST_CASE("Ideal gas conversions and equal-state flux", "[unit]") {
    Gas gas;
    State<double> w = {{1.2, 0.4, -0.1, 0.8}};
    auto u = conservative(w, gas), back = primitive(u, gas);
    auto f = hllc(w, w, gas);
    for (int k = 0; k < 4; ++k) {
        REQUIRE(back[k] == Approx(w[k]));
        REQUIRE(f.value[k] == Approx(physical_flux(w, gas)[k]));
    }
    REQUIRE(sound_speed(w, gas) == Approx(std::sqrt(1.4 * 0.8 / 1.2)));
    REQUIRE_FALSE(f.fallback);
    State<double> vacuum = {{-1, 0, 0, 1}};
    REQUIRE_FALSE(physical(vacuum, gas));
    auto rotated = hllc(w, w, gas, 0.0, 1.0);
    REQUIRE(rotated.value[0] == Approx(w[0] * w[2]));
}
TEST_CASE("Limiters preserve linear slopes and suppress extrema", "[unit]") {
    for (auto l : {Limiter::Minmod, Limiter::VanLeer, Limiter::MC}) {
        REQUIRE(limit(2.0, 2.0, l) == Approx(2));
        REQUIRE(limit(-2.0, -2.0, l) == Approx(-2));
        REQUIRE(limit(-1.0, 1.0, l) == 0);
        REQUIRE(limit(0.0, 1.0, l) == 0);
    }
    REQUIRE(limit(1.0, 3.0, Limiter::Minmod) == 1);
    REQUIRE(limit(1.0, 3.0, Limiter::VanLeer) == 1.5);
    REQUIRE(limit(1.0, 3.0, Limiter::MC) == 2);
}
TEST_CASE("HLLC contact and expansion fallback", "[unit]") {
    Gas gas;
    auto f = hllc(State<double>{{1, 0, 0, 1}}, State<double>{{0.125, 0, 0, 1}}, gas);
    REQUIRE(f.value[0] == Approx(0).margin(1e-15));
    REQUIRE(f.value[1] == Approx(1));
    auto rare = hllc(State<double>{{1, -10, 0, 0.1}}, State<double>{{1, 10, 0, 0.1}}, gas);
    REQUIRE(rare.fallback);
    for (int k = 0; k < 4; ++k)
        REQUIRE(std::isfinite(rare.value[k]));
}
TEST_CASE("Uniform flow is preserved and invalid input rejected", "[unit]") {
    Euler1D<double> s(64);
    s.initialize(std::vector<State<double>>(64, {{1, 0.4, 0.2, 1}}));
    auto before = s.state();
    s.advance(0.1);
    for (std::size_t i = 0; i < before.size(); ++i)
        REQUIRE(s.state()[i] == Approx(before[i]).margin(1e-14));
    REQUIRE_THROWS(Euler1D<double>(0));
    REQUIRE_THROWS(s.initialize(std::vector<State<double>>(64, {{-1, 0, 0, 1}})));
}
TEST_CASE("Sod exact solution and conservative mass", "[verification]") {
    constexpr int n = 400;
    Euler1D<double> s(n);
    s.sod();
    s.advance(0.2);
    reference::Riemann exact({{1, 0, 0, 1}}, {{0.125, 0, 0, 0.1}}, 1.4);
    REQUIRE(exact.pstar == Approx(0.303130178).epsilon(1e-8));
    REQUIRE(exact.ustar == Approx(0.92745262).epsilon(1e-8));
    double error[4] = {}, mass = 0;
    for (int i = 0; i < n; ++i) {
        State<double> u;
        for (int k = 0; k < 4; ++k)
            u[k] = s.state()[k * n + i];
        auto w = primitive(u, Gas{}), ref = exact.sample((i + 0.5) / n, 0.2);
        REQUIRE(physical(w, Gas{}));
        mass += w[0] / n;
        for (int k = 0; k < 4; ++k)
            error[k] += std::abs(w[k] - ref[k]) / n;
    }
    std::cout << "Sod nx=400 L1 rho=" << error[0] << " u=" << error[1] << " p=" << error[3]
              << " mass_error=" << std::abs(mass - 0.5625) << '\n';
    REQUIRE(error[0] < 0.006);
    REQUIRE(error[1] < 0.008);
    REQUIRE(error[3] < 0.005);
    REQUIRE(mass == Approx(0.5625).margin(1e-12));
    REQUIRE(s.diagnostics().reconstruction_fallbacks == 0);
}

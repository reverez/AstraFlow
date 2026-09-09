#include "astraflow/io/output.hpp"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
using namespace astraflow;
using Catch::Approx;
namespace {
nlohmann::json quantities(double flow = 2) {
    nlohmann::json e;
    for (auto name : stability_fields)
        e[name] = flow;
    e["mass_conservation_error"] = 0.0001;
    e["mass_flow_spread"] = 0.0002;
    return e;
}
} // namespace
TEST_CASE("Steady monitor requires a complete stable residual window", "[unit][convergence]") {
    ConvergenceSettings s;
    s.minimum_iterations = 40;
    s.sampling_interval = 10;
    s.window_iterations = 20;
    ConvergenceMonitor m(s, 1e-6);
    m.sample(10, {1e-2, 1e-2, 1e-2, 1e-2}, quantities());
    m.sample(20, {1e-7, 1e-7, 1e-7, 1e-7}, quantities());
    m.sample(30, {1e-7, 1e-7, 1e-7, 1e-7}, quantities());
    REQUIRE_FALSE(m.converged());
    m.sample(40, {1e-7, 1e-7, 1e-7, 1e-7}, quantities());
    REQUIRE(m.converged());
    auto r = m.report({1e-7, 1e-7, 1e-7, 1e-7});
    REQUIRE(r["residual_reduction_factors"][0].get<double>() == Approx(1e5));
    REQUIRE(r["convergence_window"]["start_iteration"] == 20);
    m.sample(50, {1e-7, 1e-7, 1e-7, 2e-6}, quantities());
    REQUIRE_FALSE(m.converged());
    REQUIRE_THROWS(m.sample(50, {0, 0, 0, 0}, quantities()));
}
TEST_CASE("Observable range deviation and fitted drift use the entire window",
          "[unit][convergence]") {
    ConvergenceSettings s;
    s.minimum_iterations = 1;
    s.sampling_interval = 10;
    s.window_iterations = 20;
    ConvergenceMonitor m(s, 1e-6);
    m.sample(10, {1e-7, 1e-7, 1e-7, 1e-7}, quantities(1));
    m.sample(20, {1e-7, 1e-7, 1e-7, 1e-7}, quantities(2));
    m.sample(30, {1e-7, 1e-7, 1e-7, 1e-7}, quantities(3));
    auto r = m.report({1e-7, 1e-7, 1e-7, 1e-7})["observable_stability"]["mass_flow"];
    REQUIRE(r["relative_range"].get<double>() == Approx(1));
    REQUIRE(r["relative_stddev"].get<double>() == Approx(std::sqrt(2.0 / 3) / 2));
    REQUIRE(r["relative_drift"].get<double>() == Approx(1));
    REQUIRE_FALSE(m.converged());
}
TEST_CASE("Conservation limits and undefined observables cannot pass", "[unit][convergence]") {
    ConvergenceSettings s;
    s.minimum_iterations = 1;
    s.sampling_interval = 10;
    s.window_iterations = 10;
    for (auto key : {"mass_conservation_error", "mass_flow_spread", "specific_impulse"}) {
        ConvergenceMonitor m(s, 1e-6);
        auto e = quantities();
        if (std::string(key) == "specific_impulse")
            e[key] = nullptr;
        else
            e[key] = 0.01;
        m.sample(10, {1e-7, 1e-7, 1e-7, 1e-7}, e);
        m.sample(20, {1e-7, 1e-7, 1e-7, 1e-7}, e);
        REQUIRE_FALSE(m.converged());
    }
    s.relative_residual_target = 0.01;
    ConvergenceMonitor m(s, 1e-6);
    m.sample(10, {1e-7, 1e-7, 1e-7, 1e-7}, quantities());
    m.sample(20, {1e-7, 1e-7, 1e-7, 1e-7}, quantities());
    REQUIRE_FALSE(m.converged());
    REQUIRE_THROWS(m.sample(30, {NAN, 0, 0, 0}, quantities()));
}
TEST_CASE("Simulation limits and user stop report nonconvergence", "[unit][convergence]") {
    Config c;
    c.backend = "cpu";
    c.precision = "double";
    c.geometry.nx = 8;
    c.geometry.nr = 4;
    c.max_iterations = 2;
    Simulation s(c);
    s.step();
    s.step();
    REQUIRE(s.termination() == "iteration_limit");
    REQUIRE(s.convergence_report()["converged"] == false);
    REQUIRE(s.convergence_report()["initial_residuals"].is_array());
    c.end_time = 1e-10;
    c.max_iterations = 10;
    Simulation timed(c);
    timed.step();
    REQUIRE(timed.termination() == "time_limit");
    Simulation stopped(c);
    stopped.stop();
    REQUIRE(stopped.termination() == "user_stop");
    REQUIRE(stopped.convergence_report()["converged"] == false);
    auto j = c.json();
    j["convergence"]["window_iterations"] = 21;
    REQUIRE_THROWS(Config::parse(j));
}

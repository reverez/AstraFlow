#include "astraflow/analysis/balance.hpp"
#include "astraflow/io/output.hpp"
#include "astraflow/numerics/finite_volume.hpp"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <fstream>
using namespace astraflow;
using Catch::Approx;
TEST_CASE("Local pseudo timestep is finite and has convective and viscous scaling",
          "[unit][steady]") {
    auto coarse = rectangular_mesh(8, 8, 1, 1, true, true);
    auto fine = rectangular_mesh(16, 16, 1, 1, true, true);
    Settings s;
    State<double> w{{1, .3, .1, 1}};
    auto step = [&](const Mesh &m) { return fv_timestep(w, m.cells[0], m.faces.data(), s); };
    auto hc = step(coarse), hf = step(fine);
    REQUIRE(hc > 0);
    REQUIRE(std::isfinite(hc));
    REQUIRE(hf == Approx(.5 * hc).epsilon(1e-13));
    s.gas.viscosity = .1;
    auto vc = step(coarse), vf = step(fine);
    REQUIRE(vc > 0);
    REQUIRE(vc < hc);
    REQUIRE(vf < hf);
    REQUIRE((1 / vf - 1 / hf) / (1 / vc - 1 / hc) == Approx(4).epsilon(1e-13));
    REQUIRE(vf / vc > .25);
    REQUIRE(vf / vc < .5);
    s.cfl *= .5;
    REQUIRE(step(coarse) == Approx(vc * .5).epsilon(1e-13));
    Nozzle n;
    n.nx = 16;
    n.nr = 8;
    auto mesh = nozzle_mesh(n);
    for (const auto &c : mesh.cells) {
        auto dt = fv_timestep(w, c, mesh.faces.data(), s);
        REQUIRE(dt > 0);
        REQUIRE(std::isfinite(dt));
    }
}
TEST_CASE("Pseudo configuration is optional and validates its own time semantics",
          "[unit][steady]") {
    Config c;
    c.backend = "cpu";
    c.precision = "double";
    REQUIRE_FALSE(c.steady.enabled);
    auto j = c.json();
    j["steady_solver"] = {{"enabled", true}, {"local_time_step", true}, {"cfl", .4}};
    auto parsed = Config::parse(j);
    REQUIRE(parsed.steady.enabled);
    REQUIRE(parsed.json() == j);
    for (double bad : {0., -.1, .5001, double(INFINITY)}) {
        auto wrong = j;
        wrong["steady_solver"]["cfl"] = bad;
        REQUIRE_THROWS(Config::parse(wrong));
    }
    auto wrong = j;
    wrong["steady_solver"]["unknown"] = true;
    REQUIRE_THROWS(Config::parse(wrong));
    wrong = j;
    wrong["runtime"]["end_time"] = 1.;
    REQUIRE_THROWS(Config::parse(wrong));
    wrong = j;
    wrong["numerics"]["profile"] = true;
    REQUIRE_THROWS(Config::parse(wrong));
    wrong = j;
    wrong["steady_solver"]["enabled"] = "yes";
    REQUIRE_THROWS(Config::parse(wrong));
}
TEST_CASE("Physical and pseudo modes use identical spatial residuals", "[unit][steady]") {
    Nozzle n;
    n.nx = 16;
    n.nr = 8;
    auto m = nozzle_mesh(n);
    Settings s;
    s.no_slip = true;
    s.gas.viscosity = .0001;
    auto w = nozzle_initial_state(m, n, s);
    SteadySettings steady;
    steady.enabled = true;
    auto physical = make_cpu_solver(m, s, true), pseudo = make_cpu_solver(m, s, true, steady);
    physical->initialize(w);
    pseudo->initialize(w);
    // A common tiny cap removes only the timestep difference; compare independently advanced
    // states.
    physical->step(1e-7);
    pseudo->step(1e-7);
    auto a = physical->state(), b = pseudo->state();
    for (std::size_t i = 0; i < a.size(); ++i)
        REQUIRE(a[i] == Approx(b[i]).margin(2e-14));
    auto exact = finite_volume_balance(m, b, s);
    for (int k = 0; k < 4; ++k)
        REQUIRE(pseudo->stats().residual[k] ==
                Approx(exact.global["rms"][k].get<double>()).margin(1e-13));
    REQUIRE(pseudo->stats().time == 0);
    REQUIRE(pseudo->stats().dt == 0);
    REQUIRE(physical->stats().time > 0);
    pseudo->step();
    REQUIRE(pseudo->stats().pseudo_dt_max > pseudo->stats().pseudo_dt_min);
    exact = finite_volume_balance(m, pseudo->state(), s);
    for (int k = 0; k < 4; ++k)
        REQUIRE(pseudo->stats().residual[k] ==
                Approx(exact.global["rms"][k].get<double>()).margin(1e-13));
}
TEST_CASE("Local pseudo RK preserves fixed points and rejects invalid states", "[unit][steady]") {
    auto m = rectangular_mesh(8, 8, 1, 1, true, true);
    Settings s;
    s.gas.viscosity = .01;
    SteadySettings steady;
    steady.enabled = true;
    auto solver = make_cpu_solver(m, s, true, steady);
    std::vector<State<double>> w(m.cells.size(), {{1, .2, .1, 1}});
    solver->initialize(w);
    auto before = solver->state();
    for (int i = 0; i < 20; ++i)
        solver->step();
    auto after = solver->state();
    for (std::size_t i = 0; i < before.size(); ++i)
        REQUIRE(after[i] == Approx(before[i]).margin(1e-14));
    for (auto r : solver->stats().residual)
        REQUIRE(r < 1e-13);
    REQUIRE(solver->stats().time == 0);
    REQUIRE_THROWS(solver->step(0));
    w[0][0] = -1;
    REQUIRE_THROWS(solver->initialize(w));
}
TEST_CASE("Pseudo metadata preserves zero physical time and monitor compatibility",
          "[unit][steady]") {
    Config c;
    c.backend = "cpu";
    c.precision = "double";
    c.geometry.nx = 8;
    c.geometry.nr = 4;
    c.steady.enabled = true;
    c.problem = "viscous_channel";
    c.settings.p0 = 1;
    c.settings.t0 = 1;
    c.settings.back_pressure = 1;
    c.max_iterations = 10;
    c.convergence.minimum_iterations = 3;
    c.convergence.sampling_interval = 1;
    c.convergence.window_iterations = 2;
    std::vector<State<double>> w(32, {{1, .2, 0, 1}});
    c.output = std::filesystem::temp_directory_path() /
               ("astraflow-c2c-test-" +
                std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    Simulation sim(c, w);
    RunOutput output(sim);
    while (!sim.finished()) {
        sim.step();
        output.record(sim);
    }
    REQUIRE(sim.termination() == "steady_converged");
    REQUIRE(sim.stats().time == 0);
    REQUIRE(sim.stats().dt == 0);
    auto report = output.finish(sim, 1);
    REQUIRE(report["time_integration"] == "pseudo_time");
    REQUIRE(report["simulated_time"] == 0);
    REQUIRE(report["physical_simulated_time"] == 0);
    REQUIRE(report["pseudo_iterations"] == sim.stats().iterations);
    REQUIRE(report["mass_conservation_error"] ==
            report["interior_inlet_exit_mass_estimate_mismatch"]);
    REQUIRE(report.contains("numerical_boundary_mass_flux_mismatch"));
    nlohmann::json metadata;
    std::ifstream in(c.output / "metadata.json");
    in >> metadata;
    REQUIRE(metadata["time_integration"] == "pseudo_time");
    REQUIRE(metadata["residual_sampling"] == "completed_pseudo_iterate");
    // The legacy mass estimate still blocks an otherwise residual-closed state.
    ConvergenceMonitor monitor(c.convergence, 1e-6);
    auto engineering = sim.sampled_engineering();
    engineering["mass_conservation_error"] = .002;
    engineering["numerical_boundary_mass_flux_mismatch"] = 0;
    for (int i = 1; i <= 4; ++i)
        monitor.sample(i, {0, 0, 0, 0}, engineering);
    REQUIRE_FALSE(monitor.converged());
    std::filesystem::remove_all(c.output); // Only this test's unique temporary artifact directory.
}
#ifdef ASTRAFLOW_HAS_CUDA
TEST_CASE("CUDA local pseudo updates and final residual match CPU FP64", "[cuda][steady]") {
    Nozzle n;
    n.nx = 16;
    n.nr = 8;
    auto m = nozzle_mesh(n);
    Settings s;
    s.no_slip = true;
    s.gas.viscosity = .0001;
    SteadySettings steady;
    steady.enabled = true;
    auto initial = nozzle_initial_state(m, n, s);
    auto cpu = make_cpu_solver(m, s, true, steady), gpu = make_cuda_solver(m, s, true, steady);
    cpu->initialize(initial);
    gpu->initialize(initial);
    for (int i = 0; i < 100; ++i) {
        cpu->step();
        gpu->step();
    }
    auto a = cpu->state(), b = gpu->state();
    for (std::size_t i = 0; i < a.size(); ++i)
        REQUIRE(std::abs(a[i] - b[i]) < 2e-12);
    auto exact = finite_volume_balance(m, b, s);
    for (int k = 0; k < 4; ++k)
        REQUIRE(gpu->stats().residual[k] ==
                Approx(exact.global["rms"][k].get<double>()).margin(2e-12));
    REQUIRE(gpu->stats().time == 0);
    REQUIRE(gpu->stats().dt == 0);
    REQUIRE(gpu->stats().pseudo_dt_max > gpu->stats().pseudo_dt_min);
    REQUIRE(gpu->stats().flux_fallbacks == 0);
    REQUIRE(gpu->stats().reconstruction_fallbacks == 0);
    auto rest = rectangular_mesh(8, 8, 1, 1, true, true);
    gpu = make_cuda_solver(rest, s, true, steady);
    initial.assign(rest.cells.size(), {{1, .2, .1, 1}});
    gpu->initialize(initial);
    a = gpu->state();
    for (int i = 0; i < 20; ++i)
        gpu->step();
    b = gpu->state();
    for (std::size_t i = 0; i < a.size(); ++i)
        REQUIRE(std::abs(a[i] - b[i]) < 1e-14);
}
#endif

#include "astraflow/core/device.hpp"
#include "astraflow/geometry/nozzle.hpp"
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <nlohmann/json.hpp>
using namespace astraflow;
struct Measurement {
    double iteration_ms, wall_ms;
    std::array<double, 7> stages{};
};
Measurement measure(Solver &solver, int iterations) {
    for (int i = 0; i < 20; ++i)
        solver.step();
    double initial = solver.stats().total_ms;
    Measurement m{};
    auto start = std::chrono::steady_clock::now();
    for (int i = 0; i < iterations; ++i) {
        solver.step();
        for (int k = 0; k < 7; ++k)
            m.stages[k] += solver.stats().stage_ms[k] / iterations;
    }
    m.wall_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start)
                    .count() /
                iterations;
    m.iteration_ms = (solver.stats().total_ms - initial) / iterations;
    return m;
}
int main(int argc, char **argv) {
    try {
        int iterations = 50, nx = 0, nr = 0, threads = std::min(8, cpu_hardware_threads());
        bool profile = false;
        std::filesystem::path output = "benchmark.json";
        for (int i = 1; i < argc; ++i) {
            std::string a = argv[i];
            if (a == "--profile")
                profile = true;
            else if (i + 1 < argc && a == "--iterations")
                iterations = std::stoi(argv[++i]);
            else if (i + 1 < argc && a == "--nx")
                nx = std::stoi(argv[++i]);
            else if (i + 1 < argc && a == "--nr")
                nr = std::stoi(argv[++i]);
            else if (i + 1 < argc && a == "--output")
                output = argv[++i];
            else if (i + 1 < argc && a == "--threads")
                threads = std::stoi(argv[++i]);
            else
                throw std::invalid_argument("Usage: astraflow_benchmark [--iterations N] [--nx N "
                                            "--nr N] [--profile] [--threads N] [--output FILE]");
        }
        if (iterations < 1 || ((nx == 0) != (nr == 0)))
            throw std::invalid_argument("Invalid benchmark dimensions/count");
        if (std::filesystem::exists(output))
            throw std::runtime_error("Benchmark output already exists");
        nlohmann::json report = {
            {"device", device_info()},
            {"precision", "FP32 CPU and GPU"},
            {"astraflow_version", version()},
            {"cpu_threads", threads},
            {"omp_wait_policy",
             std::getenv("OMP_WAIT_POLICY") ? std::getenv("OMP_WAIT_POLICY") : ""},
            {"omp_proc_bind", std::getenv("OMP_PROC_BIND") ? std::getenv("OMP_PROC_BIND") : ""},
            {"iterations", iterations},
            {"warmup_iterations", 20},
            {"profile_enabled", profile},
            {"stage_names",
             {"primitive_conversion", "CFL_reduction", "gradients", "convective_viscous_faces",
              "residual_assembly", "RK_update", "diagnostic_reductions_and_single_transfer"}},
            {"results", nlohmann::json::array()}};
        std::vector<std::pair<int, int>> grids =
            nx ? std::vector<std::pair<int, int>>{{nx, nr}}
               : std::vector<std::pair<int, int>>{{128, 32}, {256, 64}, {512, 128}, {1024, 256}};
        for (auto [x, r] : grids) {
            Nozzle g;
            g.nx = x;
            g.nr = r;
            auto mesh = nozzle_mesh(g);
            Settings s;
            s.no_slip = true;
            s.gas.viscosity = 1e-4;
            s.profile = profile;
            auto initial = nozzle_initial_state(mesh, g, s);
            // Single-thread timing is the historical reference; the parallel figure uses every
            // available OpenMP worker. Fields must agree bitwise between the two.
            int workers = threads;
            set_cpu_threads(1);
            auto cpu = make_cpu_solver(mesh, s, false);
            cpu->initialize(initial);
            auto cm = measure(*cpu, iterations);
            set_cpu_threads(workers);
            auto parallel = make_cpu_solver(mesh, s, false);
            parallel->initialize(initial);
            auto pm = measure(*parallel, iterations);
            if (parallel->state() != cpu->state())
                throw std::runtime_error("Parallel CPU backend diverged from single-thread result");
            nlohmann::json row = {{"nx", x},
                                  {"nr", r},
                                  {"cells", x * r},
                                  {"cpu_iteration_ms", cm.iteration_ms},
                                  {"cpu_wall_ms", cm.wall_ms},
                                  {"cpu_threads", workers},
                                  {"cpu_parallel_iteration_ms", pm.iteration_ms},
                                  {"cpu_parallel_speedup", cm.iteration_ms / pm.iteration_ms}};
#ifdef ASTRAFLOW_HAS_CUDA
            auto gpu = make_cuda_solver(mesh, s, false);
            gpu->initialize(initial);
            auto gm = measure(*gpu, iterations);
            row["gpu_iteration_ms"] = gm.iteration_ms;
            row["gpu_wall_ms"] = gm.wall_ms;
            row["iterations_per_second"] = 1000 / gm.iteration_ms;
            row["speedup"] = cm.iteration_ms / gm.iteration_ms;
            row["device_bytes"] = gpu->stats().device_bytes;
            row["stage_ms"] = gm.stages;
            auto a = cpu->state(), b = gpu->state();
            double max = 0, l1 = 0, l2 = 0;
            for (std::size_t i = 0; i < a.size(); ++i) {
                double d = std::abs(a[i] - b[i]);
                max = std::max(max, d);
                l1 += d / a.size();
                l2 += d * d / a.size();
            }
            row["parity_max"] = max;
            row["parity_L1"] = l1;
            row["parity_L2"] = std::sqrt(l2);
            std::cout << x << 'x' << r << " CPU=" << cm.iteration_ms << " ms CPUx" << workers << "="
                      << pm.iteration_ms << " ms GPU=" << gm.iteration_ms
                      << " ms speedup=" << cm.iteration_ms / gm.iteration_ms
                      << " memory=" << gpu->stats().device_bytes << " parity_max=" << max << '\n';
#else
            std::cout << x << 'x' << r << " CPU=" << cm.iteration_ms << " ms CPUx" << workers << "="
                      << pm.iteration_ms << " ms\n";
#endif
            report["results"].push_back(row);
        }
        if (!output.parent_path().empty())
            std::filesystem::create_directories(output.parent_path());
        std::ofstream out(output);
        out.exceptions(std::ios::failbit | std::ios::badbit);
        out << report.dump(2) << '\n';
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "Benchmark: " << e.what() << '\n';
        return 1;
    }
}

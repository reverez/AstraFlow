#include "astraflow/core/device.hpp"
#include "astraflow/io/output.hpp"
#include <chrono>
#include <csignal>
#include <fstream>
namespace {
volatile std::sig_atomic_t interrupted = 0;
void interrupt(int) { interrupted = 1; }
} // namespace
#include <iostream>
#include <optional>
int main(int argc, char **argv) {
    try {
        std::filesystem::path config_path, initial_path, state_path;
        std::optional<std::string> backend, precision, output;
        std::optional<int> iterations;
        for (int i = 1; i < argc; ++i) {
            std::string arg = argv[i];
            if (arg == "--help") {
                std::cout << "AstraFlow finite-volume CFD " << astraflow::version()
                          << "\nUsage: astraflow_cli --config FILE "
                             "[--backend cpu|cuda] [--precision float|double]\n  [--output "
                             "DIRECTORY] [--max-iterations N] [--initial-state FILE]\n  "
                             "[--threads N]   CPU reference backend worker threads\n  "
                             "[--write-state FILE]   Export the final physical state (schema 1)\n  "
                             "--version       Print the version and exit\n  "
                             "--device-info   Detect GPU and "
                             "execute native architecture probe\n";
                return 0;
            }
            if (arg == "--version") {
                std::cout << "AstraFlow " << astraflow::version() << '\n';
                return 0;
            }
            if (arg == "--device-info") {
                std::cout << astraflow::device_info() << '\n';
                return 0;
            }
            if (i + 1 >= argc)
                throw std::invalid_argument("Missing value for " + arg);
            std::string value = argv[++i];
            if (arg == "--config")
                config_path = value;
            else if (arg == "--initial-state")
                initial_path = value;
            else if (arg == "--write-state")
                state_path = value;
            else if (arg == "--backend")
                backend = value;
            else if (arg == "--precision")
                precision = value;
            else if (arg == "--output")
                output = value;
            else if (arg == "--threads") {
                std::size_t used = 0;
                int n = std::stoi(value, &used);
                if (used != value.size())
                    throw std::invalid_argument("Invalid thread count");
                astraflow::set_cpu_threads(n);
            } else if (arg == "--max-iterations") {
                std::size_t used = 0;
                int n = std::stoi(value, &used);
                if (used != value.size())
                    throw std::invalid_argument("Invalid iteration count");
                iterations = n;
            } else
                throw std::invalid_argument("Unknown option: " + arg);
        }
        if (config_path.empty())
            throw std::invalid_argument("--config is required; use --help");
        auto config = astraflow::Config::read(config_path);
        if (backend)
            config.backend = *backend;
        if (precision)
            config.precision = *precision;
        if (output)
            config.output = *output;
        if (iterations)
            config.max_iterations = *iterations;
        config.validate();
        if (config.backend == "cuda")
            std::cout << astraflow::device_info() << '\n';
        std::vector<astraflow::State<double>> initial;
        nlohmann::json initialization;
        if (!state_path.empty() && std::filesystem::exists(state_path))
            throw std::runtime_error("State export file already exists: " + state_path.string());
        if (!initial_path.empty()) {
            std::ifstream in(initial_path);
            if (!in)
                throw std::runtime_error("Cannot open initial-state file");
            in >> initialization;
            initial = astraflow::import_state(initialization, config);
        }
        astraflow::Simulation simulation(config, initial);
        astraflow::RunOutput writer(simulation);
        if (!initial_path.empty()) {
            initialization.erase("physical_conservative");
            initialization["source_file"] = std::filesystem::absolute(initial_path).string();
            std::ofstream meta(writer.directory() / "initialization.json");
            meta.exceptions(std::ios::failbit | std::ios::badbit);
            meta << initialization.dump(2) << '\n';
        }
        auto start = std::chrono::steady_clock::now();
        std::signal(SIGINT, interrupt);
        std::signal(SIGTERM, interrupt);
        try {
            while (!simulation.finished() && !interrupted) {
                simulation.step();
                if (simulation.stats().iterations % config.output_interval == 0 ||
                    simulation.sampled_iteration() == simulation.stats().iterations)
                    writer.record(simulation);
            }
        } catch (const std::exception &) {
            writer.finish(simulation, std::chrono::duration<double, std::milli>(
                                          std::chrono::steady_clock::now() - start)
                                          .count());
            throw;
        }
        if (interrupted)
            simulation.stop();
        double ms =
            std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start)
                .count();
        auto result = writer.finish(simulation, ms);
        if (!state_path.empty() && simulation.termination() != "invalid_state") {
            std::ofstream out(state_path);
            out.exceptions(std::ios::failbit | std::ios::badbit);
            out << astraflow::export_state(simulation).dump() << '\n';
        }
        std::cout << "AstraFlow " << astraflow::version() << "\nBackend: " << config.backend << " ("
                  << config.precision << ")\nCPU threads: " << simulation.stats().threads
                  << "\nTime integration: " << simulation.time_integration()
                  << "\nGrid: " << simulation.mesh().nx << 'x' << simulation.mesh().nr
                  << "\nIterations: " << result["iterations"]
                  << "\nPhysical simulated time: " << result["simulated_time"]
                  << "\nWall time (ms): " << ms << "\nTermination: " << simulation.termination()
                  << "\nConverged: " << result["converged"]
                  << "\nMass flow: " << result["mass_flow"]
                  << "\nExit Mach: " << result["exit_mach"]
                  << "\nExit pressure: " << result["exit_pressure"]
                  << "\nEstimated thrust: " << result["estimated_thrust"]
                  << "\nSpecific impulse: " << result["specific_impulse"]
                  << "\nInterior inlet/exit mass estimate mismatch: "
                  << result["interior_inlet_exit_mass_estimate_mismatch"]
                  << "\nNumerical boundary mass flux mismatch: "
                  << result["numerical_boundary_mass_flux_mismatch"]
                  << "\nFinal nondimensional residuals: " << result["residuals_nondimensional"]
                  << "\nOutput: " << writer.directory() << '\n';
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "AstraFlow: " << e.what() << '\n';
        return 1;
    }
}

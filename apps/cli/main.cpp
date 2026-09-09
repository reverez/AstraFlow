#include "astraflow/core/device.hpp"
#include "astraflow/io/output.hpp"
#include <chrono>
#include <csignal>
namespace {
volatile std::sig_atomic_t interrupted = 0;
void interrupt(int) { interrupted = 1; }
} // namespace
#include <iostream>
#include <optional>
int main(int argc, char **argv) {
    try {
        std::filesystem::path config_path;
        std::optional<std::string> backend, precision, output;
        std::optional<int> iterations;
        for (int i = 1; i < argc; ++i) {
            std::string arg = argv[i];
            if (arg == "--help") {
                std::cout << "AstraFlow finite-volume CFD\nUsage: astraflow_cli --config FILE "
                             "[--backend cpu|cuda] [--precision float|double]\n  [--output "
                             "DIRECTORY] [--max-iterations N]\n  --device-info   Detect GPU and "
                             "execute native architecture probe\n";
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
            else if (arg == "--backend")
                backend = value;
            else if (arg == "--precision")
                precision = value;
            else if (arg == "--output")
                output = value;
            else if (arg == "--max-iterations") {
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
        astraflow::Simulation simulation(config);
        astraflow::RunOutput writer(simulation);
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
        std::cout << "Backend: " << config.backend << " (" << config.precision
                  << ")\nGrid: " << simulation.mesh().nx << 'x' << simulation.mesh().nr
                  << "\nIterations: " << result["iterations"]
                  << "\nSimulated time: " << result["simulated_time"] << "\nWall time (ms): " << ms
                  << "\nTermination: " << simulation.termination()
                  << "\nConverged: " << result["converged"]
                  << "\nMass flow: " << result["mass_flow"]
                  << "\nExit Mach: " << result["exit_mach"]
                  << "\nExit pressure: " << result["exit_pressure"]
                  << "\nEstimated thrust: " << result["estimated_thrust"]
                  << "\nSpecific impulse: " << result["specific_impulse"]
                  << "\nMass conservation error: " << result["mass_conservation_error"]
                  << "\nFinal nondimensional residuals: " << result["residuals_nondimensional"]
                  << "\nOutput: " << writer.directory() << '\n';
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "AstraFlow: " << e.what() << '\n';
        return 1;
    }
}

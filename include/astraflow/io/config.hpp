#pragma once
#include "astraflow/analysis/convergence.hpp"
#include "astraflow/geometry/nozzle.hpp"
#include <filesystem>
#include <nlohmann/json.hpp>
namespace astraflow {
struct Config {
    std::string problem = "rocket_nozzle", backend = "cuda", precision = "float";
    Nozzle geometry;
    Settings settings;
    ConvergenceSettings convergence;
    double length = 1, height = 1, end_time = 0, residual_tolerance = 1e-6;
    int max_iterations = 2000, output_interval = 20;
    std::filesystem::path output;
    void validate() const;
    nlohmann::json json() const;
    static Config parse(const nlohmann::json &);
    static Config read(const std::filesystem::path &);
};
} // namespace astraflow

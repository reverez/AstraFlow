#pragma once
#include "astraflow/io/config.hpp"
namespace astraflow {
struct Scales {
    double length = 1, density = 1, velocity = 1, temperature = 1, pressure = 1, time = 1;
};
class Simulation {
  public:
    explicit Simulation(Config);
    void step();
    bool finished() const;
    std::string termination() const;
    std::vector<double> state() const;
    const Mesh &mesh() const { return physical_mesh_; }
    const Config &config() const { return config_; }
    const Scales &scales() const { return scales_; }
    StepStats stats() const;
    const StepStats &normalized_stats() const { return solver_->stats(); }

  private:
    Config config_;
    Scales scales_;
    Mesh physical_mesh_;
    std::unique_ptr<Solver> solver_;
};
} // namespace astraflow

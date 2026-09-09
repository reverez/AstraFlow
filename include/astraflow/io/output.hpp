#pragma once
#include "astraflow/analysis/engineering.hpp"
#include <fstream>
namespace astraflow {
void write_vtk(const std::filesystem::path &, const Mesh &, const std::vector<double> &, Gas);
class RunOutput {
  public:
    explicit RunOutput(const Simulation &);
    void record(const Simulation &);
    nlohmann::json finish(const Simulation &, double wall_ms);
    const std::filesystem::path &directory() const { return directory_; }

  private:
    std::filesystem::path directory_;
    std::ofstream csv_, engineering_csv_;
    int last_iteration_ = -1, last_engineering_iteration_ = -1;
};
} // namespace astraflow

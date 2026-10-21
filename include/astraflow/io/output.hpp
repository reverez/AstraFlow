#pragma once
#include "astraflow/analysis/engineering.hpp"
#include <chrono>
#include <fstream>
namespace astraflow {
void write_vtk(const std::filesystem::path &, const Mesh &, const std::vector<double> &, Gas);
// Schema-1 physical conservative cell means (row-major j*nx+i) with the configuration sections
// that must match on import. Import validates size, physical state and configuration identity.
nlohmann::json export_state(const Simulation &);
std::vector<State<double>> import_state(const nlohmann::json &, const Config &);
class RunOutput {
  public:
    explicit RunOutput(const Simulation &);
    void record(const Simulation &);
    nlohmann::json finish(const Simulation &, double wall_ms);
    const std::filesystem::path &directory() const { return directory_; }

  private:
    std::filesystem::path directory_;
    std::ofstream csv_, engineering_csv_, balance_history_;
    std::chrono::steady_clock::time_point started_;
    int last_iteration_ = -1, last_engineering_iteration_ = -1, last_balance_iteration_ = -1;
};
} // namespace astraflow

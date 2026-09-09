#pragma once
#include <array>
#include <deque>
#include <nlohmann/json.hpp>
namespace astraflow {
struct ConvergenceSettings {
    int minimum_iterations = 5000, sampling_interval = 20, window_iterations = 2000;
    // Zero disables the additional relative-reduction requirement, never the absolute target.
    double relative_residual_target = 0, observable_tolerance = 5e-4;
    double mass_mismatch_tolerance = 1e-3, station_spread_tolerance = 2e-3;
    void validate() const;
};
inline constexpr std::array<const char *, 6> stability_fields = {
    "mass_flow",   "estimated_thrust", "specific_impulse",
    "throat_mach", "exit_mach",        "exit_pressure"};
class ConvergenceMonitor {
  public:
    ConvergenceMonitor(ConvergenceSettings, double absolute_target);
    void initial(const std::array<double, 4> &);
    void sample(int iteration, const std::array<double, 4> &, const nlohmann::json &engineering);
    bool converged() const { return converged_; }
    nlohmann::json report(const std::array<double, 4> &final_residuals) const;

  private:
    struct Sample {
        int iteration;
        std::array<double, 4> residual;
        std::array<double, 8> values;
    };
    ConvergenceSettings settings_;
    double target_;
    std::array<double, 4> initial_{};
    bool initialized_ = false, converged_ = false;
    std::deque<Sample> samples_;
    nlohmann::json stability_ = nlohmann::json::object();
};
} // namespace astraflow

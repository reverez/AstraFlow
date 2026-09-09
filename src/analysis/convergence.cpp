#include "astraflow/analysis/convergence.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace astraflow {
void ConvergenceSettings::validate() const {
    if (minimum_iterations < 1 || sampling_interval < 1 || window_iterations < sampling_interval ||
        window_iterations % sampling_interval != 0)
        throw std::invalid_argument(
            "Convergence window must be a positive multiple of sampling interval");
    for (double x : {relative_residual_target, observable_tolerance, mass_mismatch_tolerance,
                     station_spread_tolerance})
        if (!std::isfinite(x) || x < 0 || x >= 1)
            throw std::invalid_argument("Invalid convergence tolerance");
}
ConvergenceMonitor::ConvergenceMonitor(ConvergenceSettings s, double target)
    : settings_(s), target_(target) {
    s.validate();
    if (!std::isfinite(target) || target < 0)
        throw std::invalid_argument("Invalid absolute residual target");
}
void ConvergenceMonitor::initial(const std::array<double, 4> &r) {
    for (double x : r)
        if (!std::isfinite(x) || x < 0)
            throw std::runtime_error("Nonfinite or negative convergence residual");
    if (!initialized_) {
        initial_ = r;
        initialized_ = true;
    }
}
void ConvergenceMonitor::sample(int iteration, const std::array<double, 4> &r,
                                const nlohmann::json &e) {
    initial(r);
    if (iteration < 1 || (!samples_.empty() && iteration <= samples_.back().iteration))
        throw std::invalid_argument("Convergence samples must have increasing positive iterations");
    Sample sample{iteration, r, {}};
    for (int k = 0; k < 6; ++k)
        sample.values[k] =
            e.at(stability_fields[k]).is_number() ? e.at(stability_fields[k]).get<double>() : NAN;
    sample.values[6] = e.at("mass_conservation_error").get<double>();
    sample.values[7] = e.at("mass_flow_spread").get<double>();
    // Undefined Isp (zero mass flow) cannot establish rocket steady convergence.
    samples_.push_back(sample);
    int cutoff = iteration - settings_.window_iterations;
    while (samples_.size() > 1 && samples_[1].iteration <= cutoff)
        samples_.pop_front();
    bool ready = samples_.size() > 1 && samples_.front().iteration <= cutoff;
    bool passed = ready && iteration >= settings_.minimum_iterations && target_ > 0;
    for (const auto &s : samples_) {
        for (int k = 0; k < 4; ++k)
            passed &= s.residual[k] <= target_ &&
                      (settings_.relative_residual_target == 0 ||
                       s.residual[k] <= initial_[k] * settings_.relative_residual_target);
        passed &= std::isfinite(s.values[6]) && s.values[6] <= settings_.mass_mismatch_tolerance;
        passed &= std::isfinite(s.values[7]) && s.values[7] <= settings_.station_spread_tolerance;
    }
    stability_ = nlohmann::json::object();
    for (int k = 0; k < 8; ++k) {
        double mean = 0, xmin = samples_.front().values[k], xmax = xmin, xmean = 0;
        bool finite = true;
        for (const auto &s : samples_) {
            finite &= std::isfinite(s.values[k]);
            mean += s.values[k] / samples_.size();
            xmean += double(s.iteration - samples_.front().iteration) / samples_.size();
            xmin = std::min(xmin, s.values[k]);
            xmax = std::max(xmax, s.values[k]);
        }
        double variance = 0, covariance = 0, xx = 0;
        for (const auto &s : samples_) {
            double dy = s.values[k] - mean, dx = s.iteration - samples_.front().iteration - xmean;
            variance += dy * dy / samples_.size();
            covariance += dx * dy;
            xx += dx * dx;
        }
        double scale = std::max(std::abs(mean), 1e-30);
        double range = (xmax - xmin) / scale;
        const char *name = k < 6    ? stability_fields[k]
                           : k == 6 ? "mass_conservation_error"
                                    : "mass_flow_spread";
        stability_[name] = {
            {"mean", mean},
            {"minimum", xmin},
            {"maximum", xmax},
            {"relative_range", finite ? nlohmann::json(range) : nlohmann::json(nullptr)},
            {"relative_stddev",
             finite ? nlohmann::json(std::sqrt(variance) / scale) : nlohmann::json(nullptr)},
            {"relative_drift",
             finite ? nlohmann::json(xx > 0 ? covariance / xx *
                                                  (iteration - samples_.front().iteration) / scale
                                            : 0)
                    : nlohmann::json(nullptr)}};
        if (k < 6)
            passed &= finite && range <= settings_.observable_tolerance;
    }
    converged_ = passed;
}
nlohmann::json ConvergenceMonitor::report(const std::array<double, 4> &final) const {
    nlohmann::json factors = nlohmann::json::array();
    for (int k = 0; k < 4; ++k)
        factors.push_back(initialized_ && final[k] > 0 ? nlohmann::json(initial_[k] / final[k])
                                                       : nlohmann::json(nullptr));
    return {
        {"converged", converged_},
        {"initial_residuals", initialized_ ? nlohmann::json(initial_) : nlohmann::json(nullptr)},
        {"final_residuals", final},
        {"residual_reduction_factors", factors},
        {"convergence_window",
         {{"required_iterations", settings_.window_iterations},
          {"sampling_interval", settings_.sampling_interval},
          {"samples", samples_.size()},
          {"start_iteration", samples_.empty() ? 0 : samples_.front().iteration},
          {"end_iteration", samples_.empty() ? 0 : samples_.back().iteration},
          {"minimum_iterations", settings_.minimum_iterations},
          {"absolute_residual_target", target_},
          {"relative_residual_target", settings_.relative_residual_target},
          {"observable_relative_range_target", settings_.observable_tolerance},
          {"mass_mismatch_target", settings_.mass_mismatch_tolerance},
          {"station_spread_target", settings_.station_spread_tolerance}}},
        {"observable_stability", stability_}};
}
} // namespace astraflow

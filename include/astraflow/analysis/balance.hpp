#pragma once
#include "astraflow/core/solver.hpp"
#include <nlohmann/json.hpp>
namespace astraflow {
struct Balance {
    std::vector<State<double>> primitive, sx, sr, gx, gr, convective, transport, source, total;
    std::vector<int> reconstruction_faces, fallback_faces;
    nlohmann::json global;
};
// Snapshot RHS, not the RK-stage statistic. Coordinates/state/settings must share units.
Balance finite_volume_balance(const Mesh &, const std::vector<double> &, Settings);
std::vector<State<double>> conservative_prolong(const Mesh &, const Mesh &,
                                                const std::vector<double> &, Gas,
                                                nlohmann::json &audit);
} // namespace astraflow

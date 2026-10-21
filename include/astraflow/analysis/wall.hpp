#pragma once
#include "astraflow/core/solver.hpp"
#include <nlohmann/json.hpp>
namespace astraflow {
// One wall face: position, viscous traction on the fluid (per unit area, physical units), heat
// flux leaving the fluid into the wall, and the wall-adjacent cell pressure/temperature.
struct WallSample {
    double x, r, length, area;
    double traction_x, traction_r, heat_flux_into_wall, pressure, temperature;
    bool upper;
};
// Evaluates the production viscous flux function at every wall face of a physical-unit snapshot.
// Mesh, state and settings must share units; axisymmetric areas are per radian.
std::vector<WallSample> wall_distribution(const Mesh &, const std::vector<double> &, Settings);
// Integrated wall loads: axial viscous force exerted on the wall, heat transfer rate into the
// wall, and peak values. Full-revolution factor 2*pi is applied for axisymmetric meshes.
nlohmann::json wall_summary(const Mesh &, const std::vector<WallSample> &);
} // namespace astraflow

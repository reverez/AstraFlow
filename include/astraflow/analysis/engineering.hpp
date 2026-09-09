#pragma once
#include "astraflow/core/simulation.hpp"
namespace astraflow {
nlohmann::json engineering(const Mesh &, const std::vector<double> &, Gas, double ambient);
}

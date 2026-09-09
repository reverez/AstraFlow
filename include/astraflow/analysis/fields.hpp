#pragma once
#include "astraflow/core/solver.hpp"
namespace astraflow {
enum class Field {
    Pressure,
    Density,
    Temperature,
    Mach,
    AxialVelocity,
    RadialVelocity,
    Speed,
    TotalEnergy,
    Vorticity
};
const char *field_name(Field);
std::vector<double> scalar_field(const Mesh &, const std::vector<double> &, Gas, Field);
} // namespace astraflow

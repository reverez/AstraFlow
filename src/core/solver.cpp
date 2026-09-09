#include "astraflow/core/solver.hpp"
#include <stdexcept>
namespace astraflow {
void validate_settings(const Settings &s) {
    if (!(s.cfl > 0 && s.cfl <= 0.5) || !std::isfinite(s.cfl) || !(s.gas.gamma > 1) ||
        !std::isfinite(s.gas.gamma) || !(s.gas.gas_constant > 0) ||
        !std::isfinite(s.gas.gas_constant) || s.gas.viscosity < 0 ||
        !std::isfinite(s.gas.viscosity) || !(s.gas.prandtl > 0))
        throw std::invalid_argument("Invalid gas or numerical settings");
}
} // namespace astraflow

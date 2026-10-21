#include "astraflow/core/solver.hpp"
#include <stdexcept>
namespace astraflow {
void SteadySettings::validate() const {
    if (!(cfl > 0 && cfl <= 0.5) || !std::isfinite(cfl))
        throw std::invalid_argument("Steady CFL must be finite and in (0, 0.5]");
}
void validate_settings(const Settings &s) {
    for (double x : {s.p0, s.t0, s.back_pressure})
        if (!(x > 0) || !std::isfinite(x))
            throw std::invalid_argument("Invalid reservoir or outlet state");
    if (s.wall_temperature < 0 || !std::isfinite(s.wall_temperature) ||
        !std::isfinite(s.upper_wall_speed))
        throw std::invalid_argument("Invalid wall state");
    for (double f : {s.gas.rho_floor, s.gas.p_floor, s.gas.temperature_floor, s.gas.prandtl})
        if (!(f > 0) || !std::isfinite(f))
            throw std::invalid_argument("Invalid thermodynamic floor or Prandtl number");
    if (!(s.cfl > 0 && s.cfl <= 0.5) || !std::isfinite(s.cfl) || !(s.gas.gamma > 1) ||
        !std::isfinite(s.gas.gamma) || !(s.gas.gas_constant > 0) ||
        !std::isfinite(s.gas.gas_constant) || s.gas.viscosity < 0 ||
        !std::isfinite(s.gas.viscosity) || !(s.gas.prandtl > 0))
        throw std::invalid_argument("Invalid gas or numerical settings");
    if (s.gas.viscosity_model != ViscosityModel::Constant) {
        if (!(s.gas.reference_temperature > 0) || !std::isfinite(s.gas.reference_temperature) ||
            s.gas.sutherland_temperature < 0 || !std::isfinite(s.gas.sutherland_temperature) ||
            !std::isfinite(s.gas.viscosity_exponent) || s.gas.viscosity_exponent < 0)
            throw std::invalid_argument("Invalid temperature-dependent viscosity parameters");
    }
}
} // namespace astraflow

#include "astraflow/analysis/engineering.hpp"
#include <algorithm>
#include <limits>
#include <numbers>
#include <stdexcept>
namespace astraflow {
nlohmann::json engineering(const Mesh &m, const std::vector<double> &state, Gas gas,
                           double ambient) {
    int n = int(m.cells.size());
    if (state.size() != std::size_t(4 * n))
        throw std::invalid_argument("Engineering state size mismatch");
    double factor = m.axisymmetric ? 2 * std::numbers::pi : 1;
    auto at = [&](int i) {
        State<double> u;
        for (int k = 0; k < 4; ++k)
            u[k] = state[k * n + i];
        auto w = primitive(u, gas);
        if (!physical(w, gas))
            throw std::runtime_error("Invalid analysis state");
        return w;
    };
    std::vector<double> station_mass(m.nx), station_mach(m.nx), areas(m.nx);
    for (int i = 0; i < m.nx; ++i)
        for (int j = 0; j < m.nr; ++j) {
            int index = j * m.nx + i;
            auto w = at(index);
            double a = factor * m.cells[index].volume / m.cells[index].dx;
            station_mass[i] += w[0] * w[1] * a;
            station_mach[i] += w[1] / sound_speed(w, gas) * a;
            areas[i] += a;
        }
    int throat = 0;
    double minimum = std::numeric_limits<double>::max();
    for (int i = 0; i < m.nx; ++i) {
        station_mach[i] /= areas[i];
        if (areas[i] < minimum) {
            minimum = areas[i];
            throat = i;
        }
    }
    double mdot = 0, inlet = 0, thrust = 0, pressure = 0, temperature = 0, mach = 0, velocity = 0,
           exit_area = 0, chamber = 0, inlet_area = 0;
    for (int j = 0; j < m.nr; ++j) {
        int ei = j * m.nx + m.nx - 1, ii = j * m.nx;
        auto e = at(ei), in = at(ii);
        double a = factor * m.faces[m.cells[ei].faces[1]].area,
               ai = factor * m.faces[m.cells[ii].faces[0]].area, dm = e[0] * e[1] * a;
        mdot += dm;
        inlet += in[0] * in[1] * ai;
        thrust += (e[0] * e[1] * e[1] + e[3] - ambient) * a;
        pressure += e[3] * a;
        temperature += e[3] / (e[0] * gas.gas_constant) * a;
        mach += std::hypot(e[1], e[2]) / sound_speed(e, gas) * a;
        velocity += e[1] * dm;
        exit_area += a;
        chamber += in[3] * ai;
        inlet_area += ai;
    }
    auto [lo, hi] = std::minmax_element(station_mass.begin(), station_mass.end());
    double scale = std::max(std::abs(*lo), std::abs(*hi));
    nlohmann::json out = {
        {"mass_flow", mdot},
        {"inlet_mass_flow", inlet},
        {"throat_mach", station_mach[throat]},
        {"exit_mach", mach / exit_area},
        {"exit_pressure", pressure / exit_area},
        {"exit_temperature", temperature / exit_area},
        {"exit_area", exit_area},
        {"estimated_thrust", thrust},
        {"chamber_to_exit_pressure_ratio", (chamber / inlet_area) / (pressure / exit_area)},
        {"station_mass_flow", station_mass},
        {"station_axial_mach", station_mach},
        {"mass_flow_spread", scale > 0 ? (*hi - *lo) / scale : 0},
        {"mass_conservation_error", scale > 0 ? std::abs(mdot - inlet) / scale : 0},
        {"ideal_gas_estimate", true}};
    out["exit_axial_velocity"] =
        std::abs(mdot) > 1e-30 ? nlohmann::json(velocity / mdot) : nlohmann::json(nullptr);
    out["specific_impulse"] =
        mdot > 1e-30 ? nlohmann::json(thrust / (mdot * 9.80665)) : nlohmann::json(nullptr);
    return out;
}
} // namespace astraflow

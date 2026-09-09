#include "astraflow/analysis/fields.hpp"
#include <stdexcept>
namespace astraflow {
const char *field_name(Field f) {
    switch (f) {
    case Field::Pressure:
        return "Pressure";
    case Field::Density:
        return "Density";
    case Field::Temperature:
        return "Temperature";
    case Field::Mach:
        return "Mach number";
    case Field::AxialVelocity:
        return "Axial velocity";
    case Field::RadialVelocity:
        return "Radial velocity";
    case Field::Speed:
        return "Velocity magnitude";
    case Field::TotalEnergy:
        return "Total specific energy";
    case Field::Vorticity:
        return "Vorticity";
    }
    return "Unknown";
}
std::vector<double> scalar_field(const Mesh &mesh, const std::vector<double> &u, Gas gas,
                                 Field field) {
    int n = int(mesh.cells.size());
    if (u.size() != std::size_t(4 * n))
        throw std::invalid_argument("Field state size mismatch");
    std::vector<State<double>> w(n);
    for (int i = 0; i < n; ++i) {
        State<double> q;
        for (int k = 0; k < 4; ++k)
            q[k] = u[k * n + i];
        w[i] = primitive(q, gas);
    }
    std::vector<double> values(n);
    for (int i = 0; i < n; ++i) {
        auto p = w[i];
        switch (field) {
        case Field::Pressure:
            values[i] = p[3];
            break;
        case Field::Density:
            values[i] = p[0];
            break;
        case Field::Temperature:
            values[i] = p[3] / (p[0] * gas.gas_constant);
            break;
        case Field::Mach:
            values[i] = std::hypot(p[1], p[2]) / sound_speed(p, gas);
            break;
        case Field::AxialVelocity:
            values[i] = p[1];
            break;
        case Field::RadialVelocity:
            values[i] = p[2];
            break;
        case Field::Speed:
            values[i] = std::hypot(p[1], p[2]);
            break;
        case Field::TotalEnergy:
            values[i] = u[3 * n + i] / p[0];
            break;
        case Field::Vorticity: {
            auto c = mesh.cells[i];
            int l = c.neighbors[0], r = c.neighbors[1], b = c.neighbors[2], t = c.neighbors[3];
            double sx = (l >= 0 && r >= 0) ? 2 : 1, sr = (b >= 0 && t >= 0) ? 2 : 1;
            if (l < 0)
                l = i;
            if (r < 0)
                r = i;
            if (b < 0)
                b = i;
            if (t < 0)
                t = i;
            double dvdr = (w[t][2] - w[b][2]) / (sr * c.dr);
            values[i] = (w[r][2] - w[l][2]) / (sx * c.dx) - c.skew * dvdr -
                        (w[t][1] - w[b][1]) / (sr * c.dr);
            break;
        }
        }
    }
    return values;
}
} // namespace astraflow

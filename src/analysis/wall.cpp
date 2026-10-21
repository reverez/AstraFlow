#include "astraflow/analysis/wall.hpp"
#include "astraflow/numerics/finite_volume.hpp"
#include <cmath>
#include <numbers>
#include <stdexcept>
namespace astraflow {
std::vector<WallSample> wall_distribution(const Mesh &m, const std::vector<double> &u, Settings s) {
    m.validate();
    validate_settings(s);
    int n = int(m.cells.size());
    if (u.size() != std::size_t(4 * n))
        throw std::invalid_argument("Wall analysis state size mismatch");
    std::vector<double> w(4 * n), gx(4 * n), gr(4 * n);
    View<double> wv{w.data(), n}, gxv{gx.data(), n}, grv{gr.data(), n};
    for (int i = 0; i < n; ++i) {
        State<double> q;
        for (int k = 0; k < 4; ++k)
            q[k] = u[k * n + i];
        auto p = primitive(q, s.gas);
        if (!physical(p, s.gas))
            throw std::invalid_argument("Nonphysical wall analysis state");
        wv.set(i, p);
    }
    for (int i = 0; i < n; ++i)
        transport_gradients(i, wv, gxv, grv, m.cells.data(), m.faces.data(), s);
    std::vector<WallSample> out;
    for (int f = 0; f < int(m.faces.size()); ++f) {
        const auto &face = m.faces[f];
        if (face.boundary != Boundary::Wall)
            continue;
        int owner = face.left >= 0 ? face.left : face.right;
        // `sign` maps the stored face normal onto the outward normal of the fluid domain.
        double sign = face.left >= 0 ? 1 : -1;
        State<double> v = s.gas.viscosity > 0
                              ? viscous_flux(f, wv, gxv, grv, m.cells.data(), m.faces.data(), s)
                              : State<double>{};
        auto p = wv.get(owner);
        WallSample sample{};
        sample.x = face.x;
        sample.r = face.r;
        sample.area = face.area;
        sample.length = m.axisymmetric ? (face.r > 0 ? face.area / face.r : 0) : face.area;
        // viscous_flux returns F_v.n_face with F_v = (0, tau, tau.u + k grad T). The viscous
        // traction on the fluid is tau.n_out. The energy component also contains the work done
        // by a moving wall (tau.u_wall), which is removed so that heat_flux_into_wall is the
        // conductive flux -k grad(T).n_out alone. Only Cartesian upper walls can move.
        sample.traction_x = sign * v[1];
        sample.traction_r = sign * v[2];
        double wall_speed = (!m.axisymmetric && face.right < 0) ? s.upper_wall_speed : 0;
        sample.heat_flux_into_wall = -sign * v[3] + wall_speed * sample.traction_x;
        sample.pressure = p[3];
        sample.temperature = p[3] / (p[0] * s.gas.gas_constant);
        sample.upper = face.right < 0;
        out.push_back(sample);
    }
    return out;
}
nlohmann::json wall_summary(const Mesh &m, const std::vector<WallSample> &samples) {
    double factor = m.axisymmetric ? 2 * std::numbers::pi : 1;
    double force = 0, heat = 0, peak_shear = 0, peak_heat = 0, area = 0;
    for (const auto &w : samples) {
        force -= factor * w.area * w.traction_x; // reaction on the wall
        heat += factor * w.area * w.heat_flux_into_wall;
        area += factor * w.area;
        peak_shear = std::max(peak_shear, std::hypot(w.traction_x, w.traction_r));
        peak_heat = std::max(peak_heat, std::abs(w.heat_flux_into_wall));
    }
    return {{"wall_faces", samples.size()},
            {"wall_area", area},
            {"wall_axial_viscous_force", force},
            {"wall_heat_transfer_rate", heat},
            {"wall_peak_shear_stress", peak_shear},
            {"wall_peak_heat_flux", peak_heat},
            {"sign_convention",
             "force: viscous reaction on the wall along +x; heat: rate leaving the fluid into "
             "the wall"}};
}
} // namespace astraflow

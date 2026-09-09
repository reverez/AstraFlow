#include "astraflow/geometry/nozzle.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace astraflow {
void Nozzle::validate() const {
    for (double v : {chamber_radius, chamber_length, throat_radius, contraction_length, exit_radius,
                     expansion_length})
        if (!(v > 0) || !std::isfinite(v))
            throw std::invalid_argument("Nozzle radii and lengths must be positive and finite");
    if (chamber_radius <= throat_radius || exit_radius <= throat_radius || nx < 4 || nr < 2 ||
        nx > 1000000 / nr)
        throw std::invalid_argument("Invalid nozzle radii or grid resolution");
}
double Nozzle::radius(double x) const {
    if (x <= chamber_length)
        return chamber_radius;
    bool contraction = x < chamber_length + contraction_length;
    double t = contraction ? (x - chamber_length) / contraction_length
                           : (x - chamber_length - contraction_length) / expansion_length;
    t = std::clamp(t, 0.0, 1.0);
    double a = contraction ? chamber_radius : throat_radius,
           b = contraction ? throat_radius : exit_radius;
    return a + (b - a) * t * t * (3 - 2 * t);
}
double Nozzle::derivative(double x) const {
    if (x <= chamber_length || x >= length())
        return 0;
    bool contraction = x < chamber_length + contraction_length;
    double l = contraction ? contraction_length : expansion_length;
    double t = (x - chamber_length - (contraction ? 0 : contraction_length)) / l;
    return (contraction ? throat_radius - chamber_radius : exit_radius - throat_radius) * 6 * t *
           (1 - t) / l;
}
Mesh nozzle_mesh(const Nozzle &g) {
    g.validate();
    Mesh m = rectangular_mesh(g.nx, g.nr, g.length(), 1);
    m.axisymmetric = true;
    int nx = g.nx, nr = g.nr;
    double dx = g.length() / nx;
    for (auto &p : m.vertices)
        p.r *= g.radius(p.x);
    for (int j = 0; j < nr; ++j)
        for (int i = 0; i <= nx; ++i) {
            auto &f = m.faces[j * (nx + 1) + i];
            double R = g.radius(f.x), lo = double(j) / nr * R, hi = double(j + 1) / nr * R;
            f.r = (lo + hi) / 2;
            f.area = (hi * hi - lo * lo) / 2;
            if (i == 0)
                f.boundary = Boundary::Inlet;
            if (i == nx)
                f.boundary = Boundary::Outlet;
        }
    int offset = (nx + 1) * nr;
    for (int j = 0; j <= nr; ++j)
        for (int i = 0; i < nx; ++i) {
            auto &f = m.faces[offset + j * nx + i];
            double eta = double(j) / nr, r0 = eta * g.radius(i * dx),
                   r1 = eta * g.radius((i + 1) * dx), dr = r1 - r0, len = std::hypot(dx, dr);
            f.r = (r0 + r1) / 2;
            f.nx = -dr / len;
            f.nr = dx / len;
            f.area = f.r * len;
            if (j == 0)
                f.boundary = Boundary::Axis;
            if (j == nr)
                f.boundary = Boundary::Wall;
        }
    for (int j = 0; j < nr; ++j)
        for (int i = 0; i < nx; ++i) {
            auto &c = m.cells[j * nx + i];
            double r0 = g.radius(i * dx), r1 = g.radius((i + 1) * dx), eta0 = double(j) / nr,
                   eta1 = double(j + 1) / nr, R = (r0 + r1) / 2;
            c.r = (eta0 + eta1) / 2 * R;
            c.dr = R / nr;
            c.skew = (eta0 + eta1) / 2 * (r1 - r0) / dx;
            c.volume = (eta1 * eta1 - eta0 * eta0) * dx * (r0 * r0 + r0 * r1 + r1 * r1) / 6;
            c.planar_area = dx * (eta1 - eta0) * R;
            c.radial_source = c.planar_area / c.volume;
        }
    m.validate();
    return m;
}
double isentropic_mach(double ratio, double gamma, bool sup) {
    if (!(ratio >= 1) || !(gamma > 1))
        throw std::invalid_argument("Invalid area-Mach relation input");
    double lo = sup ? 1 : 1e-9, hi = sup ? 20 : 1;
    for (int i = 0; i < 80; ++i) {
        double mid = (lo + hi) / 2;
        double a = std::pow(2 / (gamma + 1) * (1 + (gamma - 1) / 2 * mid * mid),
                            (gamma + 1) / (2 * (gamma - 1))) /
                   mid;
        if ((a < ratio) == sup)
            lo = mid;
        else
            hi = mid;
    }
    return (lo + hi) / 2;
}
std::vector<State<double>> nozzle_initial_state(const Mesh &m, const Nozzle &g, Settings s) {
    std::vector<State<double>> w;
    for (auto c : m.cells) {
        double R = g.radius(c.x),
               mach = isentropic_mach(R * R / (g.throat_radius * g.throat_radius), s.gas.gamma,
                                      c.x > g.chamber_length + g.contraction_length),
               factor = 1 + (s.gas.gamma - 1) / 2 * mach * mach;
        double temperature = s.t0 / factor,
               p = s.p0 / std::pow(factor, s.gas.gamma / (s.gas.gamma - 1)),
               rho = p / (s.gas.gas_constant * temperature),
               u = mach * std::sqrt(s.gas.gamma * s.gas.gas_constant * temperature);
        w.push_back({{rho, u, u * c.r / R * g.derivative(c.x), p}});
    }
    return w;
}
} // namespace astraflow

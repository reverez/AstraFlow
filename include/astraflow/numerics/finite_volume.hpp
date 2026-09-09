#pragma once
#include "astraflow/core/solver.hpp"
namespace astraflow {
template <class T> AF_HD State<T> boundary_state(State<T> w, const Face &f, Settings s) {
    if (f.boundary == Boundary::Axis) {
        w[2] = -w[2];
    }
    if (f.boundary == Boundary::Wall) {
        if (s.no_slip) {
            w[1] = -w[1];
            w[2] = -w[2];
        } else {
            T normal = w[1] * T(f.nx) + w[2] * T(f.nr);
            w[1] -= T(2 * f.nx) * normal;
            w[2] -= T(2 * f.nr) * normal;
        }
    }
    if (f.boundary == Boundary::Inlet) {
        T gamma = T(s.gas.gamma), j = w[1] - T(2) / (gamma - T(1)) * sound_speed(w, s.gas);
        T lo = 0, hi = T(0.999);
        for (int i = 0; i < (sizeof(T) == 8 ? 52 : 28); ++i) {
            T m = T(0.5) * (lo + hi);
            T a = std::sqrt(gamma * T(s.gas.gas_constant * s.t0) /
                            (T(1) + (gamma - T(1)) * T(0.5) * m * m));
            if (a * (m - T(2) / (gamma - T(1))) < j)
                lo = m;
            else
                hi = m;
        }
        T mach = T(0.5) * (lo + hi), factor = T(1) + (gamma - T(1)) * T(0.5) * mach * mach;
        T temperature = T(s.t0) / factor, p = T(s.p0) / std::pow(factor, gamma / (gamma - T(1)));
        w = {{p / (T(s.gas.gas_constant) * temperature),
              mach * std::sqrt(gamma * T(s.gas.gas_constant) * temperature), 0, p}};
    }
    if (f.boundary == Boundary::Outlet && w[1] < sound_speed(w, s.gas)) {
        T olda = sound_speed(w, s.gas),
          rho = w[0] * std::pow(T(s.back_pressure) / w[3], T(1 / s.gas.gamma));
        T a = std::sqrt(T(s.gas.gamma * s.back_pressure) / rho);
        w = {{rho, w[1] + T(2 / (s.gas.gamma - 1)) * (olda - a), w[2], T(s.back_pressure)}};
    }
    return w;
}
template <class T>
AF_HD State<T> neighbor(View<T> w, const Cell &c, int d, const Face *faces, Settings s) {
    int n = c.neighbors[d];
    if (n >= 0)
        return w.get(n);
    const auto &f = faces[c.faces[d]];
    int i = f.left >= 0 ? f.left : f.right;
    return boundary_state(w.get(i), f, s);
}
template <class T>
AF_HD void slopes(int i, View<T> w, View<T> sx, View<T> sr, const Cell *cells, const Face *faces,
                  Settings s) {
    auto c = cells[i];
    auto center = w.get(i), l = neighbor(w, c, 0, faces, s), r = neighbor(w, c, 1, faces, s),
         b = neighbor(w, c, 2, faces, s), t = neighbor(w, c, 3, faces, s);
    State<T> x, y;
    for (int k = 0; k < 4; ++k) {
        x[k] = limit(center[k] - l[k], r[k] - center[k], s.limiter);
        y[k] = limit(center[k] - b[k], t[k] - center[k], s.limiter);
    }
    sx.set(i, x);
    sr.set(i, y);
}
template <class T>
AF_HD Flux<T> fv_flux(int i, View<T> w, View<T> sx, View<T> sr, const Face *faces, Settings s,
                      bool &corrected) {
    auto f = faces[i];
    View<T> slope = f.direction ? sr : sx;
    State<T> l, r;
    if (f.left >= 0) {
        l = w.get(f.left);
        auto q = l + T(0.5) * slope.get(f.left);
        if (physical(q, s.gas))
            l = q;
        else
            corrected = true;
    }
    if (f.right >= 0) {
        r = w.get(f.right);
        auto q = r - T(0.5) * slope.get(f.right);
        if (physical(q, s.gas))
            r = q;
        else
            corrected = true;
    }
    if (f.left < 0)
        l = boundary_state(r, f, s);
    if (f.right < 0)
        r = boundary_state(l, f, s);
    return hllc(l, r, s.gas, T(f.nx), T(f.nr));
}
template <class T>
AF_HD State<T> fv_residual(int i, View<T> flux, View<T> w, const Cell *cells, const Face *faces) {
    const auto &c = cells[i];
    State<T> res;
    for (int d = 0; d < 4; ++d) {
        int f = c.faces[d];
        res = res + T((d % 2 == 0 ? 1 : -1) * faces[f].area / c.volume) * flux.get(f);
    }
    res[2] += T(c.radial_source) * w.get(i)[3];
    return res;
}
template <class T> AF_HD T fv_timestep(State<T> w, const Cell &c, const Face *faces, Settings s) {
    T sigma = 0, a = sound_speed(w, s.gas);
    for (int d = 0; d < 4; ++d) {
        auto f = faces[c.faces[d]];
        sigma += (absolute(w[1] * T(f.nx) + w[2] * T(f.nr)) + a) * T(f.area);
    }
    return T(2 * s.cfl * c.volume) / sigma;
}
} // namespace astraflow

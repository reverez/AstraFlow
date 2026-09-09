#pragma once
#include "astraflow/core/solver.hpp"
namespace astraflow {
template <class T> AF_HD State<T> boundary_state(State<T> w, const Face &f, Settings s) {
    if (f.boundary == Boundary::Axis) {
        w[2] = -w[2];
    }
    if (f.boundary == Boundary::Wall) {
        if (s.no_slip) {
            w[1] = T(2 * (f.right < 0 ? s.upper_wall_speed : 0)) - w[1];
            w[2] = -w[2];
        } else {
            T normal = w[1] * T(f.nx) + w[2] * T(f.nr);
            w[1] -= T(2 * f.nx) * normal;
            w[2] -= T(2 * f.nr) * normal;
        }
    }
    if (f.boundary == Boundary::Wall && s.wall_temperature > 0) {
        T temperature = T(2 * s.wall_temperature) - w[3] / (w[0] * T(s.gas.gas_constant));
        w[0] = w[3] / (T(s.gas.gas_constant) * temperature);
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

template <class T> AF_HD State<T> transport_variables(State<T> w, Gas g) {
    w[3] = w[3] / (w[0] * T(g.gas_constant));
    return w;
}
template <class T>
AF_HD void transport_gradients(int i, View<T> w, View<T> gx, View<T> gr, const Cell *cells,
                               const Face *faces, Settings s) {
    auto c = cells[i];
    auto l = transport_variables(neighbor(w, c, 0, faces, s), s.gas),
         r = transport_variables(neighbor(w, c, 1, faces, s), s.gas),
         b = transport_variables(neighbor(w, c, 2, faces, s), s.gas),
         t = transport_variables(neighbor(w, c, 3, faces, s), s.gas);
    auto radial = T(0.5 / c.dr) * (t - b), axial = T(0.5 / c.dx) * (r - l) - T(c.skew) * radial;
    gx.set(i, axial);
    gr.set(i, radial);
}
template <class T>
AF_HD State<T> viscous_flux(int index, View<T> w, View<T> gx, View<T> gr, const Cell *cells,
                            const Face *faces, Settings s) {
    auto f = faces[index];
    int li = f.left, ri = f.right, owner = li >= 0 ? li : ri;
    auto l = w.get(li >= 0 ? li : ri), r = w.get(ri >= 0 ? ri : li);
    if (li < 0)
        l = boundary_state(r, f, s);
    if (ri < 0)
        r = boundary_state(l, f, s);
    auto x = gx.get(owner), y = gr.get(owner);
    T dx, dr;
    if (li >= 0 && ri >= 0) {
        x = T(0.5) * (gx.get(li) + gx.get(ri));
        y = T(0.5) * (gr.get(li) + gr.get(ri));
        dx = T(cells[ri].x - cells[li].x);
        dr = T(cells[ri].r - cells[li].r);
        if (f.direction == 0 && absolute(dx) > T(1.5 * cells[owner].dx)) {
            dx = T(cells[owner].dx);
            dr = T(cells[owner].skew) * dx;
        }
        if (f.direction == 1 && absolute(dr) > T(1.5 * cells[owner].dr)) {
            dx = 0;
            dr = T(cells[owner].dr);
        }
    } else {
        T sign = li < 0 ? T(-2) : T(2);
        dx = sign * T(f.x - cells[owner].x);
        dr = sign * T(f.r - cells[owner].r);
    }
    auto jump = transport_variables(r, s.gas) - transport_variables(l, s.gas);
    T distance2 = dx * dx + dr * dr;
    for (int k = 1; k < 4; ++k) {
        T correction = (jump[k] - x[k] * dx - y[k] * dr) / distance2;
        x[k] += correction * dx;
        y[k] += correction * dr;
    }
    T u = T(0.5) * (l[1] + r[1]), v = T(0.5) * (l[2] + r[2]);
    if (f.boundary == Boundary::Wall && s.wall_temperature == 0) {
        T normal = x[3] * T(f.nx) + y[3] * T(f.nr);
        x[3] -= normal * T(f.nx);
        y[3] -= normal * T(f.nr);
    }
    T hoop = cells[owner].radial_source > 0 ? (f.r > 0 ? v / T(f.r) : y[2]) : T(0);
    T mu = T(s.gas.viscosity), div = x[1] + y[2] + hoop;
    T xx = mu * (T(2) * x[1] - T(2.0 / 3) * div), rr = mu * (T(2) * y[2] - T(2.0 / 3) * div),
      xr = mu * (y[1] + x[2]);
    T tx = xx * T(f.nx) + xr * T(f.nr), tr = xr * T(f.nx) + rr * T(f.nr);
    T k = mu * T(s.gas.gamma * s.gas.gas_constant / ((s.gas.gamma - 1) * s.gas.prandtl));
    return {{0, tx, tr, u * tx + v * tr + k * (x[3] * T(f.nx) + y[3] * T(f.nr))}};
}
template <class T>
AF_HD Flux<T> fv_flux(int i, View<T> w, View<T> sx, View<T> sr, View<T> gx, View<T> gr,
                      const Cell *cells, const Face *faces, Settings s, bool &corrected) {
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
    auto result = hllc(l, r, s.gas, T(f.nx), T(f.nr));
    if (s.gas.viscosity > 0)
        result.value = result.value - viscous_flux(i, w, gx, gr, cells, faces, s);
    return result;
}
template <class T>
AF_HD State<T> fv_residual(int i, View<T> flux, View<T> w, View<T> gx, View<T> gr,
                           const Cell *cells, const Face *faces, Settings s) {
    const auto &c = cells[i];
    State<T> res;
    for (int d = 0; d < 4; ++d) {
        int f = c.faces[d];
        res = res + T((d % 2 == 0 ? 1 : -1) * faces[f].area / c.volume) * flux.get(f);
    }
    T theta = 0;
    if (c.radial_source > 0 && s.gas.viscosity > 0) {
        T vr = w.get(i)[2] / T(c.r);
        theta = T(s.gas.viscosity) * (T(2) * vr - T(2.0 / 3) * (gx.get(i)[1] + gr.get(i)[2] + vr));
    }
    res[2] += T(c.radial_source) * (w.get(i)[3] - theta);
    return res;
}
template <class T> AF_HD T fv_timestep(State<T> w, const Cell &c, const Face *faces, Settings s) {
    T sigma = 0, a = sound_speed(w, s.gas);
    for (int d = 0; d < 4; ++d) {
        auto f = faces[c.faces[d]];
        sigma += (absolute(w[1] * T(f.nx) + w[2] * T(f.nr)) + a) * T(f.area);
    }
    T convection = T(2 * s.cfl * c.volume) / sigma;
    if (s.gas.viscosity <= 0)
        return convection;
    T diffusivity = T(s.gas.viscosity) * larger(T(4.0 / 3), T(s.gas.gamma / s.gas.prandtl)) / w[0];
    T diffusion = T(0.5 * s.cfl) /
                  (diffusivity * T(1 / (c.dx * c.dx) + (2 + c.skew * c.skew) / (c.dr * c.dr)));
    return T(1) / (T(1) / convection + T(1) / diffusion);
}
} // namespace astraflow

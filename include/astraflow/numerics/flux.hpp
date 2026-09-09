#pragma once
#include "astraflow/physics/state.hpp"
namespace astraflow {
enum class Limiter { Minmod, VanLeer, MC };
template <class T> AF_HD T minmod(T a, T b) {
    return a * b <= 0 ? T(0) : (a > 0 ? smaller(a, b) : larger(a, b));
}
template <class T> AF_HD T limit(T a, T b, Limiter limiter) {
    if (limiter == Limiter::Minmod)
        return minmod(a, b);
    if (limiter == Limiter::VanLeer)
        return a * b <= 0 ? T(0) : T(2) * a * b / (a + b);
    return minmod(T(0.5) * (a + b), minmod(T(2) * a, T(2) * b));
}
template <class T> AF_HD State<T> physical_flux(State<T> w, Gas gas) {
    auto u = conservative(w, gas);
    return {{u[1], u[1] * w[1] + w[3], u[1] * w[2], (u[3] + w[3]) * w[1]}};
}
template <class T> struct Flux {
    State<T> value;
    bool fallback = false;
};
// Rotate to the face frame; the tangential momentum is retained across contact waves.
template <class T>
AF_HD Flux<T> hllc(State<T> left, State<T> right, Gas gas, T nx = T(1), T ny = T(0)) {
    State<T> l = {{left[0], left[1] * nx + left[2] * ny, -left[1] * ny + left[2] * nx, left[3]}};
    State<T> r = {
        {right[0], right[1] * nx + right[2] * ny, -right[1] * ny + right[2] * nx, right[3]}};
    T al = sound_speed(l, gas), ar = sound_speed(r, gas);
    T sl = smaller(l[1] - al, r[1] - ar), sr = larger(l[1] + al, r[1] + ar);
    auto ul = conservative(l, gas), ur = conservative(r, gas), fl = physical_flux(l, gas),
         fr = physical_flux(r, gas);
    State<T> f;
    bool fallback = false;
    if (sl >= 0)
        f = fl;
    else if (sr <= 0)
        f = fr;
    else {
        T dl = l[0] * (sl - l[1]), dr = r[0] * (sr - r[1]);
        T sm = (r[3] - l[3] + dl * l[1] - dr * r[1]) / (dl - dr);
        T ps = l[3] + dl * (sm - l[1]);
        State<T> ls = {{dl / (sl - sm), dl / (sl - sm) * sm, dl / (sl - sm) * l[2],
                        ((sl - l[1]) * ul[3] - l[3] * l[1] + ps * sm) / (sl - sm)}};
        State<T> rs = {{dr / (sr - sm), dr / (sr - sm) * sm, dr / (sr - sm) * r[2],
                        ((sr - r[1]) * ur[3] - r[3] * r[1] + ps * sm) / (sr - sm)}};
        fallback = !std::isfinite(sm) || ps <= T(gas.p_floor) ||
                   !physical(primitive(ls, gas), gas) || !physical(primitive(rs, gas), gas);
        if (fallback)
            f = (T(1) / (sr - sl)) * (sr * fl - sl * fr + (sl * sr) * (ur - ul));
        else
            f = sm >= 0 ? fl + sl * (ls - ul) : fr + sr * (rs - ur);
    }
    T fn = f[1], ft = f[2];
    f[1] = nx * fn - ny * ft;
    f[2] = ny * fn + nx * ft;
    return {f, fallback};
}
} // namespace astraflow

#pragma once
#include "astraflow/numerics/flux.hpp"
namespace astraflow {
template <class T> AF_HD State<T> slope_1d(View<T> w, int i, Limiter limiter) {
    State<T> s{};
    if (i == 0 || i == w.size - 1)
        return s;
    auto l = w.get(i - 1), c = w.get(i), r = w.get(i + 1);
    for (int k = 0; k < 4; ++k)
        s[k] = limit(c[k] - l[k], r[k] - c[k], limiter);
    return s;
}
template <class T>
AF_HD Flux<T> face_1d(View<T> w, int f, Gas gas, Limiter limiter, bool &corrected) {
    int li = f == 0 ? 0 : f - 1, ri = f == w.size ? w.size - 1 : f;
    auto l = w.get(li), r = w.get(ri);
    if (f > 0 && f < w.size) {
        auto a = l + T(0.5) * slope_1d(w, li, limiter), b = r - T(0.5) * slope_1d(w, ri, limiter);
        if (physical(a, gas) && physical(b, gas)) {
            l = a;
            r = b;
        } else
            corrected = true;
    }
    return hllc(l, r, gas);
}
} // namespace astraflow

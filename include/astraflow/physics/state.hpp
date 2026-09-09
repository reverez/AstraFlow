#pragma once
#include <cmath>
#ifdef __CUDACC__
#define AF_HD __host__ __device__
#else
#define AF_HD
#endif
namespace astraflow {
template <class T> AF_HD T smaller(T a, T b) { return a < b ? a : b; }
template <class T> AF_HD T larger(T a, T b) { return a > b ? a : b; }
template <class T> AF_HD T absolute(T a) { return a < 0 ? -a : a; }
template <class T> struct State {
    T q[4]{};
    AF_HD T &operator[](int k) { return q[k]; }
    AF_HD T operator[](int k) const { return q[k]; }
};
template <class T> AF_HD State<T> operator+(State<T> a, State<T> b) {
    for (int k = 0; k < 4; ++k)
        a[k] += b[k];
    return a;
}
template <class T> AF_HD State<T> operator-(State<T> a, State<T> b) {
    for (int k = 0; k < 4; ++k)
        a[k] -= b[k];
    return a;
}
template <class T> AF_HD State<T> operator*(T s, State<T> a) {
    for (int k = 0; k < 4; ++k)
        a[k] *= s;
    return a;
}
struct Gas {
    double gamma = 1.4, gas_constant = 1.0, prandtl = 0.72, viscosity = 0;
    double rho_floor = 1e-10, p_floor = 1e-10, temperature_floor = 1e-10;
};
template <class T> AF_HD State<T> conservative(State<T> w, Gas g) {
    return {{w[0], w[0] * w[1], w[0] * w[2],
             w[3] / T(g.gamma - 1) + T(0.5) * w[0] * (w[1] * w[1] + w[2] * w[2])}};
}
template <class T> AF_HD State<T> primitive(State<T> u, Gas g) {
    T x = u[1] / u[0], y = u[2] / u[0];
    return {{u[0], x, y, T(g.gamma - 1) * (u[3] - T(0.5) * u[0] * (x * x + y * y))}};
}
template <class T> AF_HD bool physical(State<T> w, Gas g) {
    for (int k = 0; k < 4; ++k)
        if (!std::isfinite(w[k]))
            return false;
    return w[0] > T(g.rho_floor) && w[3] > T(g.p_floor) &&
           w[3] / (w[0] * T(g.gas_constant)) > T(g.temperature_floor);
}
template <class T> AF_HD T sound_speed(State<T> w, Gas g) { return sqrt(T(g.gamma) * w[3] / w[0]); }
template <class T> struct View {
    T *data = nullptr;
    int size = 0;
    AF_HD State<T> get(int i) const {
        State<T> s;
        for (int k = 0; k < 4; ++k)
            s[k] = data[k * size + i];
        return s;
    }
    AF_HD void set(int i, State<T> s) const {
        for (int k = 0; k < 4; ++k)
            data[k * size + i] = s[k];
    }
};
} // namespace astraflow

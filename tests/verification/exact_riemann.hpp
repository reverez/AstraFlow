#pragma once
#include "astraflow/physics/state.hpp"
#include <stdexcept>
namespace reference {
// Exact self-similar ideal-gas Riemann solution: Bracketed solve of the two wave curves.
struct Riemann {
    astraflow::State<double> l, r;
    double g, pstar, ustar;
    double wave(double p, astraflow::State<double> w) const {
        double a = std::sqrt(g * w[3] / w[0]);
        if (p > w[3])
            return (p - w[3]) * std::sqrt(2 / ((g + 1) * w[0]) / (p + (g - 1) / (g + 1) * w[3]));
        return 2 * a / (g - 1) * (std::pow(p / w[3], (g - 1) / (2 * g)) - 1);
    }
    Riemann(astraflow::State<double> left, astraflow::State<double> right, double gamma)
        : l(left), r(right), g(gamma) {
        double lo = 1e-14, hi = std::max(l[3], r[3]);
        auto f = [&](double p) { return wave(p, l) + wave(p, r) + r[1] - l[1]; };
        while (f(hi) < 0)
            hi *= 2;
        for (int i = 0; i < 100; ++i) {
            double p = (lo + hi) / 2;
            if (f(p) > 0)
                hi = p;
            else
                lo = p;
        }
        pstar = (lo + hi) / 2;
        ustar = (l[1] + r[1] + wave(pstar, r) - wave(pstar, l)) / 2;
    }
    astraflow::State<double> sample(double x, double t) const {
        double s = (x - 0.5) / t;
        bool left = s <= ustar;
        auto w = left ? l : r;
        double sign = left ? -1 : 1;
        double a = std::sqrt(g * w[3] / w[0]), ratio = pstar / w[3], beta = (g - 1) / (g + 1);
        if (pstar > w[3]) {
            double speed =
                w[1] + sign * a * std::sqrt((g + 1) / (2 * g) * ratio + (g - 1) / (2 * g));
            if (sign * (s - speed) > 0)
                return w;
            return {{w[0] * (ratio + beta) / (beta * ratio + 1), ustar, 0, pstar}};
        }
        double head = w[1] + sign * a, astar = a * std::pow(ratio, (g - 1) / (2 * g)),
               tail = ustar + sign * astar;
        if (sign * (s - head) >= 0)
            return w;
        if (sign * (s - tail) <= 0)
            return {{w[0] * std::pow(ratio, 1 / g), ustar, 0, pstar}};
        double u = 2 / (g + 1) * (-sign * a + (g - 1) / 2 * w[1] + s);
        double af = 2 / (g + 1) * (a + sign * (g - 1) / 2 * (s - w[1]));
        return {
            {w[0] * std::pow(af / a, 2 / (g - 1)), u, 0, w[3] * std::pow(af / a, 2 * g / (g - 1))}};
    }
};
} // namespace reference

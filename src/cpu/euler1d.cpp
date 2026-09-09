#include "astraflow/core/euler1d.hpp"
#include <chrono>
#include <limits>
#include <stdexcept>
#include <string>
namespace astraflow {
template <class T>
Euler1D<T>::Euler1D(int cells, Gas gas, T cfl, Limiter limiter)
    : cells_(cells), gas_(gas), cfl_(cfl), limiter_(limiter) {
    if (cells < 4 || !(cfl > 0 && cfl <= T(0.5)) || !(gas.gamma > 1) || !(gas.gas_constant > 0))
        throw std::invalid_argument("Euler1D requires nx>=4, 0<CFL<=0.5 and valid ideal gas");
    u_.resize(4 * cells);
    w_.resize(4 * cells);
    stage_.resize(4 * cells);
    residual_.resize(4 * cells);
    fluxes_.resize(4 * (cells + 1));
}
template <class T> void Euler1D<T>::initialize(const std::vector<State<T>> &w) {
    if (w.size() != static_cast<std::size_t>(cells_))
        throw std::invalid_argument("Initial state size mismatch");
    View<T> u{u_.data(), cells_};
    for (int i = 0; i < cells_; ++i) {
        if (!physical(w[i], gas_))
            throw std::invalid_argument("Nonphysical initial state");
        u.set(i, conservative(w[i], gas_));
    }
    diagnostics_ = {};
}
template <class T> void Euler1D<T>::sod() {
    std::vector<State<T>> w(cells_);
    for (int i = 0; i < cells_; ++i)
        w[i] = i < cells_ / 2 ? State<T>{{1, 0, 0, 1}} : State<T>{{T(0.125), 0, 0, T(0.1)}};
    initialize(w);
}
template <class T> void Euler1D<T>::validate(const std::vector<T> &state) const {
    for (int i = 0; i < cells_; ++i) {
        State<T> u;
        for (int k = 0; k < 4; ++k)
            u[k] = state[k * cells_ + i];
        if (!physical(primitive(u, gas_), gas_))
            throw std::runtime_error("Invalid Euler state at cell " + std::to_string(i) +
                                     ", iteration " + std::to_string(diagnostics_.iterations));
    }
}
template <class T> T Euler1D<T>::timestep() const {
    T dt = std::numeric_limits<T>::max();
    for (int i = 0; i < cells_; ++i) {
        State<T> u;
        for (int k = 0; k < 4; ++k)
            u[k] = u_[k * cells_ + i];
        auto w = primitive(u, gas_);
        dt = smaller(dt, cfl_ / (T(cells_) * (absolute(w[1]) + sound_speed(w, gas_))));
    }
    return dt;
}
template <class T> void Euler1D<T>::rhs(const std::vector<T> &state) {
    View<T> w{w_.data(), cells_}, flux{fluxes_.data(), cells_ + 1}, res{residual_.data(), cells_};
    for (int i = 0; i < cells_; ++i) {
        State<T> u;
        for (int k = 0; k < 4; ++k)
            u[k] = state[k * cells_ + i];
        w.set(i, primitive(u, gas_));
    }
    for (int f = 0; f <= cells_; ++f) {
        bool corrected = false;
        auto value = face_1d(w, f, gas_, limiter_, corrected);
        flux.set(f, value.value);
        diagnostics_.flux_fallbacks += value.fallback;
        diagnostics_.reconstruction_fallbacks += corrected;
    }
    for (int i = 0; i < cells_; ++i)
        res.set(i, T(cells_) * (flux.get(i) - flux.get(i + 1)));
}
template <class T> void Euler1D<T>::advance(T end_time) {
    if (!std::isfinite(end_time) || end_time < T(diagnostics_.time))
        throw std::invalid_argument("Invalid end time");
    auto start = std::chrono::steady_clock::now();
    validate(u_);
    while (diagnostics_.time < double(end_time)) {
        T dt = smaller(timestep(), T(double(end_time) - diagnostics_.time));
        if (!(dt > 0) || !std::isfinite(dt))
            throw std::runtime_error("Invalid Euler time step");
        rhs(u_);
        for (std::size_t k = 0; k < u_.size(); ++k)
            stage_[k] = u_[k] + dt * residual_[k];
        validate(stage_);
        rhs(stage_);
        for (std::size_t k = 0; k < u_.size(); ++k)
            u_[k] = T(0.5) * (u_[k] + stage_[k] + dt * residual_[k]);
        validate(u_);
        diagnostics_.time += double(dt);
        ++diagnostics_.iterations;
    }
    diagnostics_.elapsed_ms +=
        std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
}
template class Euler1D<float>;
template class Euler1D<double>;
} // namespace astraflow

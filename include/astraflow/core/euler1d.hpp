#pragma once
#include "astraflow/numerics/euler1d_ops.hpp"
#include <cstdint>
#include <vector>
namespace astraflow {
struct Diagnostics {
    std::uint64_t flux_fallbacks = 0, reconstruction_fallbacks = 0;
    int iterations = 0;
    double time = 0, elapsed_ms = 0;
};
template <class T> class Euler1D {
  public:
    Euler1D(int cells, Gas gas = {}, T cfl = T(0.4), Limiter limiter = Limiter::MC);
    void initialize(const std::vector<State<T>> &primitives);
    void sod();
    void advance(T end_time);
    const std::vector<T> &state() const { return u_; }
    const Diagnostics &diagnostics() const { return diagnostics_; }

  private:
    void rhs(const std::vector<T> &state);
    T timestep() const;
    void validate(const std::vector<T> &state) const;
    int cells_;
    Gas gas_;
    T cfl_;
    Limiter limiter_;
    std::vector<T> u_, w_, stage_, residual_, fluxes_;
    Diagnostics diagnostics_;
};
extern template class Euler1D<float>;
extern template class Euler1D<double>;
} // namespace astraflow

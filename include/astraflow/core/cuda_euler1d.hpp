#pragma once
#include "astraflow/core/euler1d.hpp"
#include <memory>
namespace astraflow {
template <class T> class CudaEuler1D {
  public:
    CudaEuler1D(int cells, Gas gas = {}, T cfl = T(0.4), Limiter limiter = Limiter::MC);
    ~CudaEuler1D();
    CudaEuler1D(const CudaEuler1D &) = delete;
    CudaEuler1D &operator=(const CudaEuler1D &) = delete;
    void initialize(const std::vector<State<T>> &primitives);
    void sod();
    void advance(T end_time);
    std::vector<T> state() const;
    const Diagnostics &diagnostics() const;

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
extern template class CudaEuler1D<float>;
extern template class CudaEuler1D<double>;
} // namespace astraflow

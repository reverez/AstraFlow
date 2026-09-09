#pragma once
#include <cuda_runtime.h>
#include <stdexcept>
#include <string>
#include <vector>
namespace astraflow::cuda {
inline void check(cudaError_t e, const char *expression, int line) {
    if (e != cudaSuccess)
        throw std::runtime_error(std::string(expression) + " at line " + std::to_string(line) +
                                 ": " + cudaGetErrorString(e));
}
#define AF_CUDA(call) ::astraflow::cuda::check((call), #call, __LINE__)
template <class T> class Buffer {
  public:
    explicit Buffer(std::size_t count) : count_(count) {
        AF_CUDA(cudaMalloc(&data_, count * sizeof(T)));
    }
    ~Buffer() {
        if (data_)
            cudaFree(data_);
    }
    Buffer(const Buffer &) = delete;
    Buffer &operator=(const Buffer &) = delete;
    T *data() const { return data_; }
    std::size_t bytes() const { return count_ * sizeof(T); }
    void zero() { AF_CUDA(cudaMemset(data_, 0, bytes())); }
    void upload(const T *values) {
        AF_CUDA(cudaMemcpy(data_, values, bytes(), cudaMemcpyHostToDevice));
    }
    std::vector<T> download() const {
        std::vector<T> v(count_);
        AF_CUDA(cudaMemcpy(v.data(), data_, bytes(), cudaMemcpyDeviceToHost));
        return v;
    }

  private:
    T *data_ = nullptr;
    std::size_t count_;
};
class Event {
  public:
    Event() { AF_CUDA(cudaEventCreate(&event_)); }
    ~Event() { cudaEventDestroy(event_); }
    Event(const Event &) = delete;
    Event &operator=(const Event &) = delete;
    void record() { AF_CUDA(cudaEventRecord(event_)); }
    float since(const Event &start) {
        AF_CUDA(cudaEventSynchronize(event_));
        float ms = 0;
        AF_CUDA(cudaEventElapsedTime(&ms, start.event_, event_));
        return ms;
    }

  private:
    cudaEvent_t event_{};
};
} // namespace astraflow::cuda

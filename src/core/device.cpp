#include "astraflow/core/device.hpp"
#ifndef ASTRAFLOW_HAS_CUDA
namespace astraflow {
std::string device_info() { return "CPU reference build; CUDA disabled"; }
} // namespace astraflow
#endif

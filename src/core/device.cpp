#include "astraflow/core/device.hpp"
#ifndef ASTRAFLOW_VERSION
#define ASTRAFLOW_VERSION "unknown"
#endif
namespace astraflow {
const char *version() { return ASTRAFLOW_VERSION; }
#ifndef ASTRAFLOW_HAS_CUDA
std::string device_info() { return "CPU reference build; CUDA disabled"; }
#endif
} // namespace astraflow

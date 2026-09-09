#include "astraflow/core/device.hpp"
#include <iostream>
#include <stdexcept>
#include <string_view>
int main(int argc, char **argv) {
    try {
        if (argc == 2 && std::string_view(argv[1]) == "--device-info") {
            std::cout << astraflow::device_info() << '\n';
            return 0;
        }
        if (argc == 1 || (argc == 2 && std::string_view(argv[1]) == "--help")) {
            std::cout << "AstraFlow: Milestone 0 toolchain scaffold\n"
                         "Usage: astraflow_cli [--help | --device-info]\n"
                         "Simulation commands are not implemented yet.\n";
            return 0;
        }
        throw std::invalid_argument("Unsupported argument; use --help");
    } catch (const std::exception &e) {
        std::cerr << "AstraFlow: " << e.what() << '\n';
        return 1;
    }
}

#pragma once
#include "astraflow/core/solver.hpp"
#include "astraflow/geometry/mesh.hpp"
namespace astraflow {
struct Nozzle {
    double chamber_radius = 0.2, chamber_length = 1, throat_radius = 0.15, contraction_length = 2,
           exit_radius = 0.225, expansion_length = 4;
    int nx = 96, nr = 12;
    void validate() const;
    double length() const { return chamber_length + contraction_length + expansion_length; }
    double radius(double x) const;
    double derivative(double x) const;
};
Mesh nozzle_mesh(const Nozzle &);
double isentropic_mach(double area_ratio, double gamma, bool supersonic);
std::vector<State<double>> nozzle_initial_state(const Mesh &, const Nozzle &, Settings);
} // namespace astraflow

#pragma once
#include <vector>
namespace astraflow {
struct Point {
    double x = 0, r = 0;
};
enum class Boundary { Interior, Extrapolate, Axis, Wall, Inlet, Outlet };
struct Face {
    int left = -1, right = -1, direction = 0;
    double x = 0, r = 0, nx = 0, nr = 0, area = 0;
    Boundary boundary = Boundary::Interior;
};
struct Cell {
    int faces[4]{}, neighbors[4]{};
    double x = 0, r = 0, volume = 0, planar_area = 0, dx = 0, dr = 0, skew = 0, radial_source = 0;
};
struct Mesh {
    int nx = 0, nr = 0;
    bool axisymmetric = false;
    std::vector<Point> vertices;
    std::vector<Cell> cells;
    std::vector<Face> faces;
    void validate() const;
};
Mesh rectangular_mesh(int nx, int nr, double length = 1, double height = 1, bool periodic_x = false,
                      bool periodic_r = false, bool walls = false);
} // namespace astraflow

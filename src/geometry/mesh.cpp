#include "astraflow/geometry/mesh.hpp"
#include <cmath>
#include <stdexcept>
namespace astraflow {
void Mesh::validate() const {
    if (nx < 2 || nr < 2 || cells.size() != std::size_t(nx * nr) ||
        vertices.size() != std::size_t((nx + 1) * (nr + 1)))
        throw std::invalid_argument("Invalid mesh dimensions");
    for (auto p : vertices)
        if (!std::isfinite(p.x) || !std::isfinite(p.r))
            throw std::invalid_argument("Nonfinite mesh vertex");
    for (const auto &c : cells) {
        if (!(c.volume > 0) || !std::isfinite(c.volume) || !(c.dx > 0) || !(c.dr > 0))
            throw std::invalid_argument("Invalid cell metric");
        for (int f : c.faces)
            if (f < 0 || f >= int(faces.size()))
                throw std::invalid_argument("Invalid face index");
    }
    for (const auto &f : faces) {
        if (f.area < 0 || !std::isfinite(f.area) || std::abs(f.nx * f.nx + f.nr * f.nr - 1) > 1e-12)
            throw std::invalid_argument("Invalid face metric");
        if (f.left < 0 && f.right < 0)
            throw std::invalid_argument("Face has no adjacent cell");
        if (f.left >= int(cells.size()) || f.right >= int(cells.size()))
            throw std::invalid_argument("Invalid face adjacency");
    }
}
Mesh rectangular_mesh(int nx, int nr, double length, double height, bool px, bool pr, bool walls) {
    if (nx < 2 || nr < 2 || !(length > 0) || !(height > 0) || !std::isfinite(length) ||
        !std::isfinite(height) || nx > 1000000 / nr)
        throw std::invalid_argument("Invalid rectangular mesh parameters");
    Mesh m;
    m.nx = nx;
    m.nr = nr;
    double dx = length / nx, dr = height / nr;
    for (int j = 0; j <= nr; ++j)
        for (int i = 0; i <= nx; ++i)
            m.vertices.push_back({i * dx, j * dr});
    auto cell = [&](int i, int j) { return j * nx + i; };
    auto xf = [&](int i, int j) { return j * (nx + 1) + i; };
    auto rf = [&](int i, int j) { return (nx + 1) * nr + j * nx + i; };
    for (int j = 0; j < nr; ++j)
        for (int i = 0; i <= nx; ++i) {
            Face f;
            f.left = i ? cell(i - 1, j) : (px ? cell(nx - 1, j) : -1);
            f.right = i < nx ? cell(i, j) : (px ? cell(0, j) : -1);
            f.x = i * dx;
            f.r = (j + 0.5) * dr;
            f.nx = 1;
            f.area = dr;
            f.boundary = (f.left < 0 || f.right < 0) ? Boundary::Extrapolate : Boundary::Interior;
            m.faces.push_back(f);
        }
    for (int j = 0; j <= nr; ++j)
        for (int i = 0; i < nx; ++i) {
            Face f;
            f.left = j ? cell(i, j - 1) : (pr ? cell(i, nr - 1) : -1);
            f.right = j < nr ? cell(i, j) : (pr ? cell(i, 0) : -1);
            f.x = (i + 0.5) * dx;
            f.r = j * dr;
            f.nr = 1;
            f.area = dx;
            f.direction = 1;
            f.boundary = (f.left < 0 || f.right < 0)
                             ? (walls ? Boundary::Wall : Boundary::Extrapolate)
                             : Boundary::Interior;
            m.faces.push_back(f);
        }
    for (int j = 0; j < nr; ++j)
        for (int i = 0; i < nx; ++i) {
            Cell c;
            c.x = (i + 0.5) * dx;
            c.r = (j + 0.5) * dr;
            c.dx = dx;
            c.dr = dr;
            c.volume = dx * dr;
            c.planar_area = c.volume;
            c.faces[0] = xf(i, j);
            c.faces[1] = xf(i + 1, j);
            c.faces[2] = rf(i, j);
            c.faces[3] = rf(i, j + 1);
            for (int d = 0; d < 4; ++d) {
                const auto &f = m.faces[c.faces[d]];
                c.neighbors[d] = (d % 2 == 0) ? f.left : f.right;
            }
            m.cells.push_back(c);
        }
    m.validate();
    return m;
}
} // namespace astraflow

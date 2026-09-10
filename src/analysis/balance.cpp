#include "astraflow/analysis/balance.hpp"
#include "astraflow/numerics/finite_volume.hpp"
#include <algorithm>
#include <stdexcept>
namespace astraflow {
Balance finite_volume_balance(const Mesh &m, const std::vector<double> &u, Settings s) {
    m.validate();
    validate_settings(s);
    int n = int(m.cells.size()), nf = int(m.faces.size());
    if (u.size() != std::size_t(4 * n))
        throw std::invalid_argument("Balance state size mismatch");
    Balance b;
    for (auto *v : {&b.primitive, &b.sx, &b.sr, &b.gx, &b.gr, &b.convective, &b.transport,
                    &b.source, &b.total})
        v->resize(n);
    b.reconstruction_faces.resize(n);
    b.fallback_faces.resize(n);
    std::vector<double> w(4 * n), sx(4 * n), sr(4 * n), gx(4 * n), gr(4 * n), fc(4 * nf),
        fv(4 * nf), full(4 * nf);
    View<double> wv{w.data(), n}, sxv{sx.data(), n}, srv{sr.data(), n}, gxv{gx.data(), n},
        grv{gr.data(), n}, cv{fc.data(), nf}, vv{fv.data(), nf}, all{full.data(), nf};
    for (int i = 0; i < n; ++i) {
        State<double> q;
        for (int k = 0; k < 4; ++k)
            q[k] = u[k * n + i];
        auto p = primitive(q, s.gas);
        if (!physical(p, s.gas))
            throw std::invalid_argument("Nonphysical balance state");
        wv.set(i, p);
    }
    for (int i = 0; i < n; ++i) {
        slopes(i, wv, sxv, srv, m.cells.data(), m.faces.data(), s);
        transport_gradients(i, wv, gxv, grv, m.cells.data(), m.faces.data(), s);
    }
    std::array<State<double>, 4> boundary_c{},
        boundary_v{}; // inlet,outlet,wall,axis; outward signed.
    for (int f = 0; f < nf; ++f) {
        bool corrected = false;
        Settings inv = s;
        inv.gas.viscosity = 0;
        auto c = fv_flux(f, wv, sxv, srv, gxv, grv, m.cells.data(), m.faces.data(), inv, corrected);
        auto v = s.gas.viscosity > 0
                     ? viscous_flux(f, wv, gxv, grv, m.cells.data(), m.faces.data(), s)
                     : State<double>{};
        cv.set(f, c.value);
        vv.set(f, v);
        all.set(f, c.value - v);
        auto face = m.faces[f];
        for (int i : {face.left, face.right})
            if (i >= 0) {
                b.reconstruction_faces[i] += corrected;
                b.fallback_faces[i] += c.fallback;
            }
        if (face.left >= 0 && face.right >= 0)
            continue;
        int group = face.boundary == Boundary::Inlet    ? 0
                    : face.boundary == Boundary::Outlet ? 1
                    : face.boundary == Boundary::Wall   ? 2
                                                        : 3;
        double outward = face.left < 0 ? -1 : 1;
        boundary_c[group] = boundary_c[group] + (outward * face.area) * c.value;
        boundary_v[group] = boundary_v[group] + (outward * face.area) * v;
    }
    State<double> integral{}, conv{}, visc{}, source{}, pressure_source{}, hoop_source{}, rms{},
        defect{};
    double local_defect = 0;
    for (int i = 0; i < n; ++i) {
        auto cell = m.cells[i];
        b.primitive[i] = wv.get(i);
        b.sx[i] = sxv.get(i);
        b.sr[i] = srv.get(i);
        b.gx[i] = gxv.get(i);
        b.gr[i] = grv.get(i);
        for (int d = 0; d < 4; ++d) {
            int f = cell.faces[d];
            double sign = (d % 2 == 0 ? 1 : -1) * m.faces[f].area / cell.volume;
            b.convective[i] = b.convective[i] + sign * cv.get(f);
            b.transport[i] = b.transport[i] - sign * vv.get(f);
        }
        double theta = 0;
        if (cell.radial_source > 0 && s.gas.viscosity > 0) {
            double vr = wv.get(i)[2] / cell.r;
            theta = s.gas.viscosity * (2 * vr - (2.0 / 3) * (gxv.get(i)[1] + grv.get(i)[2] + vr));
        }
        b.source[i][2] = cell.radial_source * (wv.get(i)[3] - theta);
        b.total[i] = fv_residual(i, all, wv, gxv, grv, m.cells.data(), m.faces.data(), s);
        integral = integral + cell.volume * b.total[i];
        conv = conv + cell.volume * b.convective[i];
        visc = visc + cell.volume * b.transport[i];
        source = source + cell.volume * b.source[i];
        pressure_source[2] += cell.volume * cell.radial_source * wv.get(i)[3];
        hoop_source[2] -= cell.volume * cell.radial_source * theta;
        for (int k = 0; k < 4; ++k) {
            rms[k] += b.total[i][k] * b.total[i][k] / n;
            local_defect = std::max(local_defect, std::abs(b.total[i][k] - b.convective[i][k] -
                                                           b.transport[i][k] - b.source[i][k]));
        }
    }
    auto array = [](State<double> q) { return nlohmann::json::array({q[0], q[1], q[2], q[3]}); };
    State<double> prediction = source, denom{};
    nlohmann::json boundary;
    const char *names[] = {"inlet", "outlet", "wall", "axis_or_other"};
    for (int g = 0; g < 4; ++g) {
        boundary[names[g]] = {{"outward_convective", array(boundary_c[g])},
                              {"outward_viscous_thermal", array(boundary_v[g])},
                              {"outward_total", array(boundary_c[g] - boundary_v[g])}};
        prediction = prediction - boundary_c[g] + boundary_v[g];
        for (int k = 0; k < 4; ++k)
            denom[k] += std::abs(boundary_c[g][k]) + std::abs(boundary_v[g][k]);
    }
    State<double> relative{}, scaled_rhs{};
    for (int k = 0; k < 4; ++k) {
        rms[k] = std::sqrt(rms[k]);
        defect[k] = integral[k] - prediction[k];
        denom[k] += std::abs(pressure_source[k]) + std::abs(hoop_source[k]);
        relative[k] = std::abs(defect[k]) / std::max(denom[k], 1e-30);
        scaled_rhs[k] = std::abs(integral[k]) / std::max(denom[k], 1e-30);
    }
    b.global = {{"boundary_fluxes", boundary},
                {"integrated_rhs", array(integral)},
                {"integrated_convective", array(conv)},
                {"integrated_viscous_thermal", array(visc)},
                {"integrated_source", array(source)},
                {"pressure_source", array(pressure_source)},
                {"hoop_stress_source", array(hoop_source)},
                {"boundary_source_prediction", array(prediction)},
                {"assembly_identity_error", array(defect)},
                {"relative_identity_error", array(relative)},
                {"relative_global_rhs", array(scaled_rhs)},
                {"rms", array(rms)},
                {"local_decomposition_max_error", local_defect}};
    return b;
}
std::vector<State<double>> conservative_prolong(const Mesh &coarse, const Mesh &fine,
                                                const std::vector<double> &u, Gas gas,
                                                nlohmann::json &audit) {
    coarse.validate();
    fine.validate();
    int nx = coarse.nx, nr = coarse.nr, n = nx * nr;
    if (fine.nx != 2 * nx || fine.nr != 2 * nr || u.size() != std::size_t(4 * n))
        throw std::invalid_argument("Prolongation requires nested logical r=2 grids");
    auto at = [&](int i, int j) {
        State<double> q;
        for (int k = 0; k < 4; ++k)
            q[k] = u[k * n + j * nx + i];
        return q;
    };
    std::vector<State<double>> output(fine.cells.size());
    double error = 0, scale_change = 0;
    int limited = 0;
    for (int j = 0; j < nr; ++j)
        for (int i = 0; i < nx; ++i) {
            auto q = at(i, j);
            State<double> sx{}, sr{};
            if (!physical(primitive(q, gas), gas))
                throw std::invalid_argument("Invalid parent state");
            for (int k = 0; k < 4; ++k) {
                if (i > 0 && i < nx - 1)
                    sx[k] = limit(q[k] - at(i - 1, j)[k], at(i + 1, j)[k] - q[k], Limiter::MC);
                if (j > 0 && j < nr - 1)
                    sr[k] = limit(q[k] - at(i, j - 1)[k], at(i, j + 1)[k] - q[k], Limiter::MC);
            }
            int ids[4];
            State<double> delta[4], mean{};
            double volume = 0;
            for (int b = 0; b < 2; ++b)
                for (int a = 0; a < 2; ++a) {
                    int c = 2 * b + a;
                    ids[c] = (2 * j + b) * fine.nx + 2 * i + a;
                    double v = fine.cells[ids[c]].volume;
                    delta[c] = (a ? .25 : -.25) * sx + (b ? .25 : -.25) * sr;
                    mean = mean + v * delta[c];
                    volume += v;
                }
            mean = (1 / volume) * mean;
            double ratio = coarse.cells[j * nx + i].volume / volume;
            scale_change = std::max(scale_change, std::abs(ratio - 1));
            auto base = ratio * q;
            double alpha = 1;
            State<double> children[4];
            for (;;) {
                bool positive = true;
                for (int c = 0; c < 4; ++c) {
                    children[c] = base + alpha * (delta[c] - mean);
                    positive &= physical(primitive(children[c], gas), gas);
                }
                if (positive)
                    break;
                if (alpha < 1e-12) {
                    alpha = 0;
                } else
                    alpha *= .5;
                if (alpha == 0 && !physical(primitive(base, gas), gas))
                    throw std::runtime_error("Prolongation base is nonphysical");
            }
            limited += alpha < 1;
            State<double> sum{};
            for (int c = 0; c < 4; ++c) {
                sum = sum + fine.cells[ids[c]].volume * children[c];
                output[ids[c]] = children[c];
            }
            for (int k = 0; k < 4; ++k)
                error = std::max(
                    error, std::abs(sum[k] - q[k] * coarse.cells[j * nx + i].volume) /
                               std::max(std::abs(q[k] * coarse.cells[j * nx + i].volume), 1e-30));
        }
    audit = {{"maximum_parent_relative_integral_error", error},
             {"maximum_volume_ratio_change", scale_change},
             {"positivity_limited_parents", limited},
             {"method",
              "MC-limited conservative linear prolongation with child-volume mean correction"}};
    return output;
}
} // namespace astraflow

#include "astraflow/io/config.hpp"
#include <fstream>
#include <stdexcept>
namespace astraflow {
namespace {
void keys(const nlohmann::json &j, std::initializer_list<const char *> allowed) {
    if (!j.is_object())
        throw std::invalid_argument("Expected JSON object");
    for (auto it = j.begin(); it != j.end(); ++it) {
        bool found = false;
        for (auto k : allowed)
            found |= it.key() == k;
        if (!found)
            throw std::invalid_argument("Unknown configuration key: " + it.key());
    }
}
std::string limiter_name(Limiter l) {
    return l == Limiter::MC ? "mc" : l == Limiter::Minmod ? "minmod" : "van_leer";
}
} // namespace
void Config::validate() const {
    validate_settings(settings);
    geometry.validate();
    if (problem != "sod" && problem != "rocket_nozzle" && problem != "isentropic_nozzle" &&
        problem != "viscous_channel")
        throw std::invalid_argument("Unknown problem");
    if (backend != "cpu" && backend != "cuda")
        throw std::invalid_argument("Backend must be cpu or cuda");
    if (precision != "float" && precision != "double")
        throw std::invalid_argument("Precision must be float or double");
    for (double x : {length, height, settings.p0, settings.t0, settings.back_pressure})
        if (!(x > 0) || !std::isfinite(x))
            throw std::invalid_argument(
                "Lengths and boundary pressure/temperature must be positive and finite");
    for (double x : {end_time, residual_tolerance, settings.wall_temperature})
        if (x < 0 || !std::isfinite(x))
            throw std::invalid_argument("Invalid runtime or wall temperature");
    if (!std::isfinite(settings.upper_wall_speed) || max_iterations < 1 || output_interval < 1)
        throw std::invalid_argument("Invalid iteration settings or wall velocity");
    if (problem == "rocket_nozzle" && settings.gas.viscosity > 0 && !settings.no_slip)
        throw std::invalid_argument("Viscous rocket nozzle requires no-slip wall");
}
nlohmann::json Config::json() const {
    const auto &g = settings.gas;
    const auto &n = geometry;
    return {{"problem", problem},
            {"gas",
             {{"gamma", g.gamma},
              {"gas_constant", g.gas_constant},
              {"prandtl", g.prandtl},
              {"viscosity", g.viscosity},
              {"rho_floor", g.rho_floor},
              {"p_floor", g.p_floor},
              {"temperature_floor", g.temperature_floor}}},
            {"geometry",
             {{"chamber_radius", n.chamber_radius},
              {"chamber_length", n.chamber_length},
              {"throat_radius", n.throat_radius},
              {"contraction_length", n.contraction_length},
              {"exit_radius", n.exit_radius},
              {"expansion_length", n.expansion_length},
              {"length", length},
              {"height", height}}},
            {"mesh", {{"nx", n.nx}, {"nr", n.nr}}},
            {"boundary_conditions",
             {{"P0", settings.p0},
              {"T0", settings.t0},
              {"back_pressure", settings.back_pressure},
              {"wall", settings.no_slip ? "no_slip" : "slip"},
              {"wall_temperature", settings.wall_temperature},
              {"upper_wall_speed", settings.upper_wall_speed}}},
            {"numerics",
             {{"cfl", settings.cfl},
              {"limiter", limiter_name(settings.limiter)},
              {"profile", settings.profile}}},
            {"runtime",
             {{"backend", backend},
              {"precision", precision},
              {"max_iterations", max_iterations},
              {"end_time", end_time},
              {"residual_tolerance", residual_tolerance}}},
            {"output", {{"directory", output.string()}, {"interval", output_interval}}}};
}
Config Config::parse(const nlohmann::json &j) {
    keys(j, {"problem", "gas", "geometry", "mesh", "boundary_conditions", "numerics", "runtime",
             "output"});
    Config c;
    c.problem = j.value("problem", c.problem);
    auto read = [&](const char *section, std::initializer_list<const char *> allowed) {
        auto v = j.value(section, nlohmann::json::object());
        keys(v, allowed);
        return v;
    };
    auto g = read("gas", {"gamma", "gas_constant", "prandtl", "viscosity", "rho_floor", "p_floor",
                          "temperature_floor"});
    auto &gas = c.settings.gas;
    gas.gamma = g.value("gamma", gas.gamma);
    gas.gas_constant = g.value("gas_constant", gas.gas_constant);
    gas.prandtl = g.value("prandtl", gas.prandtl);
    gas.viscosity = g.value("viscosity", gas.viscosity);
    gas.rho_floor = g.value("rho_floor", gas.rho_floor);
    gas.p_floor = g.value("p_floor", gas.p_floor);
    gas.temperature_floor = g.value("temperature_floor", gas.temperature_floor);
    auto n =
        read("geometry", {"chamber_radius", "chamber_length", "throat_radius", "contraction_length",
                          "exit_radius", "expansion_length", "length", "height"});
    auto &shape = c.geometry;
    shape.chamber_radius = n.value("chamber_radius", shape.chamber_radius);
    shape.chamber_length = n.value("chamber_length", shape.chamber_length);
    shape.throat_radius = n.value("throat_radius", shape.throat_radius);
    shape.contraction_length = n.value("contraction_length", shape.contraction_length);
    shape.exit_radius = n.value("exit_radius", shape.exit_radius);
    shape.expansion_length = n.value("expansion_length", shape.expansion_length);
    c.length = n.value("length", c.length);
    c.height = n.value("height", c.height);
    auto m = read("mesh", {"nx", "nr"});
    shape.nx = m.value("nx", shape.nx);
    shape.nr = m.value("nr", c.problem == "sod" ? 2 : shape.nr);
    auto b = read("boundary_conditions",
                  {"P0", "T0", "back_pressure", "wall", "wall_temperature", "upper_wall_speed"});
    c.settings.p0 = b.value("P0", c.settings.p0);
    c.settings.t0 = b.value("T0", c.settings.t0);
    c.settings.back_pressure = b.value("back_pressure", c.settings.back_pressure);
    c.settings.wall_temperature = b.value("wall_temperature", 0.0);
    c.settings.upper_wall_speed = b.value("upper_wall_speed", 0.0);
    auto wall = b.value("wall", std::string("slip"));
    if (wall != "slip" && wall != "no_slip")
        throw std::invalid_argument("Unknown wall type");
    c.settings.no_slip = wall == "no_slip";
    auto nu = read("numerics", {"cfl", "limiter", "profile"});
    c.settings.profile = nu.value("profile", false);
    c.settings.cfl = nu.value("cfl", c.settings.cfl);
    auto limiter = nu.value("limiter", std::string("mc"));
    if (limiter == "mc")
        c.settings.limiter = Limiter::MC;
    else if (limiter == "minmod")
        c.settings.limiter = Limiter::Minmod;
    else if (limiter == "van_leer")
        c.settings.limiter = Limiter::VanLeer;
    else
        throw std::invalid_argument("Unknown limiter");
    auto r = read("runtime",
                  {"backend", "precision", "max_iterations", "end_time", "residual_tolerance"});
    c.backend = r.value("backend", c.backend);
    c.precision = r.value("precision", c.precision);
    c.max_iterations = r.value("max_iterations", c.max_iterations);
    c.end_time = r.value("end_time", c.end_time);
    c.residual_tolerance = r.value("residual_tolerance", c.residual_tolerance);
    auto o = read("output", {"directory", "interval"});
    c.output = o.value("directory", std::string{});
    c.output_interval = o.value("interval", c.output_interval);
    c.validate();
    return c;
}
Config Config::read(const std::filesystem::path &path) {
    std::ifstream in(path);
    if (!in)
        throw std::runtime_error("Cannot open configuration: " + path.string());
    nlohmann::json j;
    in >> j;
    return parse(j);
}
} // namespace astraflow

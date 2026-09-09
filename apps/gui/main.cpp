#include "astraflow/analysis/fields.hpp"
#include "astraflow/core/controller.hpp"
#include "astraflow/core/device.hpp"
#include <GLFW/glfw3.h>
#include <algorithm>
#include <array>
#include <fstream>
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#include <implot.h>
#include <iostream>
#include <vector>
using namespace astraflow;
namespace {
void theme() {
    ImGui::StyleColorsDark();
    auto &s = ImGui::GetStyle();
    s.WindowRounding = 8;
    s.ChildRounding = 8;
    s.FrameRounding = 5;
    s.GrabRounding = 4;
    s.WindowPadding = {16, 14};
    s.FramePadding = {8, 6};
    s.ItemSpacing = {8, 9};
    s.Colors[ImGuiCol_WindowBg] = {0.035f, 0.055f, 0.085f, 1};
    s.Colors[ImGuiCol_ChildBg] = {0.055f, 0.08f, 0.115f, 1};
    s.Colors[ImGuiCol_FrameBg] = {0.09f, 0.13f, 0.18f, 1};
    s.Colors[ImGuiCol_Button] = {0.08f, 0.26f, 0.29f, 1};
    s.Colors[ImGuiCol_ButtonHovered] = {0.10f, 0.42f, 0.43f, 1};
    s.Colors[ImGuiCol_Header] = {0.08f, 0.25f, 0.29f, 1};
    s.Colors[ImGuiCol_CheckMark] = {0.2f, 0.85f, 0.73f, 1};
    s.Colors[ImGuiCol_Text] = {0.88f, 0.93f, 0.96f, 1};
}
void capture(GLFWwindow *window, const std::filesystem::path &path) {
    int w, h;
    glfwGetFramebufferSize(window, &w, &h);
    std::vector<unsigned char> pixels(std::size_t(w) * h * 3);
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glReadBuffer(GL_BACK);
    glReadPixels(0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, pixels.data());
    if (std::filesystem::exists(path))
        throw std::runtime_error("Capture file already exists");
    if (!path.parent_path().empty())
        std::filesystem::create_directories(path.parent_path());
    std::ofstream out(path, std::ios::binary);
    out.exceptions(std::ios::failbit | std::ios::badbit);
    out << "P6\n" << w << ' ' << h << "\n255\n";
    for (int y = h - 1; y >= 0; --y)
        out.write(reinterpret_cast<const char *>(pixels.data() + std::size_t(y) * w * 3), w * 3);
}
void input(const char *label, double &value, const char *format = "%.5g") {
    ImGui::SetNextItemWidth(145);
    ImGui::InputDouble(label, &value, 0, 0, format);
}
void metric(const nlohmann::json &j, const char *key, const char *label, const char *units) {
    ImGui::TableNextRow();
    ImGui::TableNextColumn();
    ImGui::TextUnformatted(label);
    ImGui::TableNextColumn();
    if (j.contains(key) && j[key].is_number())
        ImGui::Text("%.6g %s", j[key].get<double>(), units);
    else
        ImGui::TextUnformatted("--");
}
struct ViewState {
    float zoom = 1, radial_scale = 6;
    ImVec2 pan{};
    bool mesh = false;
};
void flow_view(const Snapshot &snap, const std::vector<double> &values, ViewState &view) {
    ImVec2 origin = ImGui::GetCursorScreenPos(), size = ImGui::GetContentRegionAvail();
    size.y = std::max(size.y, 60.0f);
    ImGui::InvisibleButton("flow_canvas", size, ImGuiButtonFlags_MouseButtonMiddle);
    bool hovered = ImGui::IsItemHovered();
    if (hovered) {
        view.zoom = std::clamp(view.zoom * std::pow(1.15f, ImGui::GetIO().MouseWheel), 0.3f, 15.0f);
        if (ImGui::IsMouseDragging(ImGuiMouseButton_Middle)) {
            view.pan.x += ImGui::GetIO().MouseDelta.x;
            view.pan.y += ImGui::GetIO().MouseDelta.y;
        }
    }
    auto *draw = ImGui::GetWindowDrawList();
    draw->AddRectFilled(origin, {origin.x + size.x, origin.y + size.y}, IM_COL32(9, 18, 29, 255),
                        5);
    draw->PushClipRect(origin, {origin.x + size.x, origin.y + size.y}, true);
    const auto &m = *snap.mesh;
    double xmin = m.vertices.front().x, xmax = m.vertices[m.nx].x, rmax = 0;
    for (auto p : m.vertices)
        rmax = std::max(rmax, p.r);
    bool mirror = m.axisymmetric;
    float usable = size.x - 110;
    double scale = std::min(double(usable - 50) / (xmax - xmin),
                            double(size.y - 70) / ((mirror ? 2 : 1) * rmax * view.radial_scale)) *
                   view.zoom;
    auto project = [&](Point p, double sign) {
        return ImVec2(float(origin.x + usable / 2 + (p.x - (xmin + xmax) / 2) * scale) + view.pan.x,
                      float(origin.y + size.y / 2 -
                            (p.r * sign - (mirror ? 0 : rmax / 2)) * scale * view.radial_scale) +
                          view.pan.y);
    };
    auto [lo, hi] = std::minmax_element(values.begin(), values.end());
    double span = *hi - *lo;
    for (int j = 0; j < m.nr; ++j)
        for (int i = 0; i < m.nx; ++i) {
            int index = j * m.nx + i;
            float normalized = span > 0 ? float((values[index] - *lo) / span) : 0.5f;
            ImU32 color = ImGui::ColorConvertFloat4ToU32(
                ImPlot::SampleColormap(normalized, ImPlotColormap_Viridis));
            std::array<Point, 4> p = {
                m.vertices[j * (m.nx + 1) + i], m.vertices[j * (m.nx + 1) + i + 1],
                m.vertices[(j + 1) * (m.nx + 1) + i + 1], m.vertices[(j + 1) * (m.nx + 1) + i]};
            for (int side = 0; side < (mirror ? 2 : 1); ++side) {
                ImVec2 q[4];
                for (int k = 0; k < 4; ++k)
                    q[k] = project(p[side ? 3 - k : k], side ? -1 : 1);
                draw->AddConvexPolyFilled(q, 4, color);
                if (view.mesh)
                    draw->AddPolyline(q, 4, IM_COL32(210, 230, 240, 45), ImDrawFlags_Closed, 0.5f);
            }
        }
    for (int side = 0; side < (mirror ? 2 : 1); ++side)
        for (int i = 0; i < m.nx; ++i)
            draw->AddLine(project(m.vertices[m.nr * (m.nx + 1) + i], side ? -1 : 1),
                          project(m.vertices[m.nr * (m.nx + 1) + i + 1], side ? -1 : 1),
                          IM_COL32(214, 235, 240, 220), 1.5f);
    draw->AddLine(project({xmin, 0}, 1), project({xmax, 0}, 1), IM_COL32(130, 170, 190, 120), 1);
    float bx = origin.x + size.x - 70, top = origin.y + 35, bottom = origin.y + size.y - 35;
    for (int i = 0; i < 100; ++i) {
        float a = float(i) / 100, b = float(i + 1) / 100;
        auto color =
            ImGui::ColorConvertFloat4ToU32(ImPlot::SampleColormap(1 - a, ImPlotColormap_Viridis));
        draw->AddRectFilled({bx, top + a * (bottom - top)}, {bx + 14, top + b * (bottom - top)},
                            color);
    }
    char text[64];
    std::snprintf(text, sizeof(text), "%.4g", *hi);
    draw->AddText({bx - 4, top - 20}, IM_COL32(220, 230, 240, 255), text);
    std::snprintf(text, sizeof(text), "%.4g", *lo);
    draw->AddText({bx - 4, bottom + 6}, IM_COL32(220, 230, 240, 255), text);
    draw->AddText({origin.x + 14, origin.y + 12}, IM_COL32(150, 175, 190, 255),
                  mirror ? "MERIDIONAL FIELD / MIRRORED ABOUT AXIS" : "CARTESIAN FIELD");
    draw->AddText({origin.x + 14, origin.y + size.y - 25}, IM_COL32(150, 175, 190, 255),
                  "Cell values  |  Wheel: zoom  |  Middle drag: pan");
    draw->PopClipRect();
}
} // namespace
int main(int argc, char **argv) {
    GLFWwindow *window = nullptr;
    bool gui_started = false;
    try {
        std::filesystem::path config_path = "examples/rocket_nozzle/config.json", capture_path;
        bool smoke = false;
        int frame_limit = 0;
        for (int i = 1; i < argc; ++i) {
            std::string a = argv[i];
            if (a == "--smoke-test") {
                smoke = true;
                frame_limit = 240;
            } else if (a == "--help") {
                std::cout << "astraflow_gui [--config FILE] [--smoke-test] [--capture FILE.ppm]\n";
                return 0;
            } else if (i + 1 < argc && a == "--config")
                config_path = argv[++i];
            else if (i + 1 < argc && a == "--capture")
                capture_path = argv[++i];
            else
                throw std::invalid_argument("Unknown or incomplete GUI option");
        }
        Config active = Config::read(config_path), draft = active;
        std::string gpu = device_info(), edit_error;
        glfwSetErrorCallback(
            [](int, const char *message) { std::cerr << "GLFW: " << message << '\n'; });
        if (!glfwInit())
            throw std::runtime_error(
                "GLFW initialization failed; check WSLg DISPLAY and X11 access");
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
        window = glfwCreateWindow(1480, 980, "AstraFlow | Rocket nozzle CFD", nullptr, nullptr);
        if (!window)
            throw std::runtime_error("OpenGL window creation failed");
        glfwMakeContextCurrent(window);
        glfwSwapInterval(1);
        IMGUI_CHECKVERSION();
        ImGui::CreateContext();
        ImPlot::CreateContext();
        auto &io = ImGui::GetIO();
        io.IniFilename = nullptr;
        io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
        std::filesystem::path font = "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf";
        if (std::filesystem::exists(font))
            io.Fonts->AddFontFromFileTTF(font.c_str(), 16);
        theme();
        if (!ImGui_ImplGlfw_InitForOpenGL(window, true) || !ImGui_ImplOpenGL3_Init("#version 330"))
            throw std::runtime_error("ImGui backend initialization failed");
        gui_started = true;
        Controller controller(active);
        ViewState view;
        int selected = int(Field::Mach), previous_field = -1;
        std::shared_ptr<const Snapshot> previous;
        std::vector<double> values, hx;
        std::array<std::vector<double>, 4> hy;
        double snapshot_hz = 0, last_snapshot = glfwGetTime();
        int frame = 0, paused_count = -1, smoke_actions = 0;
        std::array<bool, 9> rendered{};
        while (!glfwWindowShouldClose(window)) {
            glfwPollEvents();
            ImGui_ImplOpenGL3_NewFrame();
            ImGui_ImplGlfw_NewFrame();
            ImGui::NewFrame();
            ++frame;
            auto snapshot = controller.snapshot();
            auto state = controller.state();
            bool running = state == RunState::Running;
            if (smoke && frame < 180)
                selected = (frame / 18) % 9;
            else if (smoke)
                selected = int(Field::Mach);
            if (snapshot && (snapshot != previous || selected != previous_field)) {
                values = scalar_field(*snapshot->mesh, snapshot->conservative, snapshot->gas,
                                      Field(selected));
                rendered[selected] = true;
                if (snapshot != previous) {
                    double now = glfwGetTime();
                    snapshot_hz = 1 / std::max(now - last_snapshot, 1e-6);
                    last_snapshot = now;
                    hx.clear();
                    for (auto &h : hy)
                        h.clear();
                    for (auto row : snapshot->history) {
                        hx.push_back(row.iteration);
                        for (int k = 0; k < 4; ++k)
                            hy[k].push_back(std::max(row.residual[k], 1e-30));
                    }
                }
                previous = snapshot;
                previous_field = selected;
            }
            ImGui::SetNextWindowPos({0, 0});
            ImGui::SetNextWindowSize(io.DisplaySize);
            ImGui::Begin("AstraFlow", nullptr,
                         ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove |
                             ImGuiWindowFlags_NoSavedSettings);
            ImGui::TextColored({0.25f, 0.87f, 0.75f, 1}, "ASTRAFLOW");
            ImGui::SameLine();
            ImGui::TextDisabled(" /  AXISYMMETRIC COMPRESSIBLE CFD");
            ImGui::SameLine(ImGui::GetWindowWidth() - 230);
            ImGui::TextColored(running ? ImVec4(0.3f, 0.9f, 0.65f, 1)
                                       : ImVec4(0.8f, 0.85f, 0.9f, 1),
                               "%s", state_name(state));
            ImGui::Separator();
            ImGui::BeginChild("Controls", {326, 0}, ImGuiChildFlags_Borders);
            ImGui::TextUnformatted("SIMULATION");
            auto reset = [&](bool regenerate) {
                try {
                    if (regenerate) {
                        draft.validate();
                        active = draft;
                    } else
                        draft = active;
                    controller.reset(active);
                    edit_error.clear();
                } catch (const std::exception &e) {
                    edit_error = e.what();
                }
            };
            bool can_run = state == RunState::Ready || state == RunState::Paused;
            ImGui::BeginDisabled(!can_run);
            if (ImGui::Button("Run", {142, 0}) || (smoke && (frame == 25 || frame == 160))) {
                controller.run();
                ++smoke_actions;
            }
            ImGui::SameLine();
            if (ImGui::Button("Step", {142, 0}) || (smoke && frame == 100)) {
                if (snapshot)
                    paused_count = snapshot->stats.iterations;
                controller.single_step();
                ++smoke_actions;
            }
            ImGui::EndDisabled();
            if (ImGui::Button("Pause", {142, 0}) || (smoke && (frame == 85 || frame == 220))) {
                controller.pause();
                ++smoke_actions;
            }
            ImGui::SameLine();
            ImGui::BeginDisabled(running);
            if (ImGui::Button("Reset", {142, 0}) || (smoke && frame == 115)) {
                reset(false);
                ++smoke_actions;
            }
            ImGui::EndDisabled();
            if (smoke && frame == 105 &&
                (!snapshot || snapshot->stats.iterations != paused_count + 1))
                throw std::runtime_error("GUI smoke single-step invariant failed");
            if (smoke && frame == 125 && (!snapshot || snapshot->stats.iterations != 0))
                throw std::runtime_error("GUI smoke reset invariant failed");
            ImGui::Spacing();
            ImGui::BeginDisabled(running || state == RunState::Idle);
            int backend = draft.backend == "cuda" ? 1 : 0,
                precision = draft.precision == "double" ? 1 : 0;
            ImGui::SetNextItemWidth(145);
            if (ImGui::Combo("Backend", &backend, "CPU\0CUDA\0"))
                draft.backend = backend ? "cuda" : "cpu";
            ImGui::SetNextItemWidth(145);
            if (ImGui::Combo("Precision", &precision, "FP32\0FP64\0"))
                draft.precision = precision ? "double" : "float";
            ImGui::SetNextItemWidth(145);
            ImGui::InputInt("Axial cells", &draft.geometry.nx);
            ImGui::SetNextItemWidth(145);
            ImGui::InputInt("Radial cells", &draft.geometry.nr);
            input("CFL", draft.settings.cfl);
            ImGui::SetNextItemWidth(145);
            ImGui::InputInt("Max iterations", &draft.max_iterations);
            input("Residual target", draft.residual_tolerance);
            if (ImGui::CollapsingHeader("Gas and chamber", ImGuiTreeNodeFlags_DefaultOpen)) {
                input("Gamma", draft.settings.gas.gamma);
                input("R [J/kg/K]", draft.settings.gas.gas_constant);
                input("Viscosity [Pa s]", draft.settings.gas.viscosity);
                input("Prandtl", draft.settings.gas.prandtl);
                input("P0 [Pa]", draft.settings.p0);
                input("T0 [K]", draft.settings.t0);
                input("Back pressure [Pa]", draft.settings.back_pressure);
            }
            if (ImGui::CollapsingHeader("Nozzle geometry [m]", ImGuiTreeNodeFlags_DefaultOpen)) {
                input("Chamber radius", draft.geometry.chamber_radius);
                input("Chamber length", draft.geometry.chamber_length);
                input("Throat radius", draft.geometry.throat_radius);
                input("Contraction length", draft.geometry.contraction_length);
                input("Exit radius", draft.geometry.exit_radius);
                input("Expansion length", draft.geometry.expansion_length);
            }
            if (ImGui::Button("Regenerate Mesh", {-1, 0}) || (smoke && frame == 135)) {
                if (smoke)
                    draft.geometry.nx += 8;
                reset(true);
                ++smoke_actions;
            }
            ImGui::EndDisabled();
            ImGui::TextDisabled("Edits apply with Regenerate Mesh.");
            if (!edit_error.empty())
                ImGui::TextWrapped("%s", edit_error.c_str());
            if (state == RunState::Failed)
                ImGui::TextWrapped("%s", controller.error().c_str());
            ImGui::EndChild();
            ImGui::SameLine();
            ImGui::BeginGroup();
            ImGui::BeginChild("Field", {0, ImGui::GetContentRegionAvail().y * 0.55f},
                              ImGuiChildFlags_Borders);
            ImGui::TextUnformatted("FLOW FIELD");
            ImGui::SameLine();
            ImGui::SetNextItemWidth(200);
            if (ImGui::BeginCombo("##field", field_name(Field(selected)))) {
                for (int f = 0; f < 9; ++f)
                    if (ImGui::Selectable(field_name(Field(f)), selected == f))
                        selected = f;
                ImGui::EndCombo();
            }
            ImGui::SameLine();
            ImGui::Checkbox("Mesh", &view.mesh);
            ImGui::SameLine();
            ImGui::SetNextItemWidth(120);
            ImGui::SliderFloat("Radial display scale", &view.radial_scale, 1, 10, "%.1fx");
            ImGui::SameLine();
            if (ImGui::SmallButton("Fit")) {
                view.zoom = 1;
                view.pan = {};
            }
            if (snapshot && !values.empty())
                flow_view(*snapshot, values, view);
            else
                ImGui::TextUnformatted("Preparing solver and mesh...");
            ImGui::EndChild();
            float available = ImGui::GetContentRegionAvail().x;
            ImGui::BeginChild("Convergence", {available * 0.56f, 0}, ImGuiChildFlags_Borders);
            ImGui::TextUnformatted("CONVERGENCE");
            ImGui::TextDisabled("RMS residuals / nondimensional solver units");
            if (ImPlot::BeginPlot("##residual", {-1, ImGui::GetContentRegionAvail().y - 115})) {
                ImPlot::SetupAxes("Iteration", "Residual");
                ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
                ImPlot::SetupAxisLimits(ImAxis_X1, 0, hx.empty() ? 100 : hx.back() + 1,
                                        ImGuiCond_Always);
                ImPlot::SetupAxisLimits(ImAxis_Y1, 1e-7, 1, ImGuiCond_Once);
                const char *labels[] = {"Mass", "Axial momentum", "Radial momentum", "Energy"};
                for (int k = 0; k < 4; ++k)
                    if (!hx.empty())
                        ImPlot::PlotLine(labels[k], hx.data(), hy[k].data(), int(hx.size()));
                ImPlot::EndPlot();
            }
            if (snapshot) {
                const auto &s = snapshot->stats;
                ImGui::Text("Iteration %d / %d    Time %.6g s", s.iterations, active.max_iterations,
                            s.time);
                ImGui::ProgressBar(float(s.iterations) / active.max_iterations, {-1, 0});
                ImGui::Text("Flux fallbacks: %llu  Reconstruction: %llu",
                            static_cast<unsigned long long>(s.flux_fallbacks),
                            static_cast<unsigned long long>(s.reconstruction_fallbacks));
            }
            ImGui::EndChild();
            ImGui::SameLine();
            ImGui::BeginChild("Engineering", {0, 0}, ImGuiChildFlags_Borders);
            ImGui::TextUnformatted("ENGINEERING ANALYSIS");
            ImGui::TextDisabled("Calorically perfect ideal-gas estimate");
            if (snapshot && ImGui::BeginTable("Metrics", 2, ImGuiTableFlags_RowBg)) {
                auto &j = snapshot->engineering;
                metric(j, "mass_flow", "Mass flow", "kg/s");
                metric(j, "throat_mach", "Throat Mach", "");
                metric(j, "exit_mach", "Exit Mach", "");
                metric(j, "exit_axial_velocity", "Exit velocity", "m/s");
                metric(j, "exit_pressure", "Exit pressure", "Pa");
                metric(j, "estimated_thrust", "Estimated thrust", "N");
                metric(j, "specific_impulse", "Specific impulse", "s");
                metric(j, "mass_conservation_error", "Mass-flow mismatch", "");
                ImGui::EndTable();
            }
            ImGui::Spacing();
            ImGui::Separator();
            ImGui::TextUnformatted("GPU / PERFORMANCE");
            ImGui::TextWrapped("%s", gpu.substr(0, gpu.find('\n')).c_str());
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("%s", gpu.c_str());
            ImGui::TextDisabled("%s", gpu.substr(gpu.find_last_of('\n') + 1).c_str());
            if (snapshot) {
                auto &s = snapshot->stats;
                ImGui::Text("Grid: %d x %d   |   %.2f MiB device", snapshot->mesh->nx,
                            snapshot->mesh->nr, double(s.device_bytes) / (1024 * 1024));
                ImGui::Text("%.3f ms/iteration  |  %.0f iterations/s", s.iteration_ms,
                            s.iteration_ms > 0 ? 1000 / s.iteration_ms : 0);
                ImGui::Text("Snapshot updates: %.1f Hz", snapshot_hz);
            }
            ImGui::EndChild();
            ImGui::EndGroup();
            ImGui::End();
            ImGui::Render();
            int width, height;
            glfwGetFramebufferSize(window, &width, &height);
            glViewport(0, 0, width, height);
            glClearColor(0.035f, 0.055f, 0.085f, 1);
            glClear(GL_COLOR_BUFFER_BIT);
            ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
            if (frame_limit && frame >= frame_limit) {
                if (!capture_path.empty())
                    capture(window, capture_path);
                glfwSwapBuffers(window);
                break;
            }
            glfwSwapBuffers(window);
            if (smoke && controller.state() == RunState::Failed)
                throw std::runtime_error(controller.error());
        }
        if (smoke) {
            if (smoke_actions != 7 ||
                !std::all_of(rendered.begin(), rendered.end(), [](bool x) { return x; }))
                throw std::runtime_error("GUI smoke did not exercise all actions/fields");
            std::cout << "GUI smoke PASS: window, nine fields, run/pause/step/reset/regenerate, "
                         "residuals and engineering dashboard\n";
        }
        ImGui_ImplOpenGL3_Shutdown();
        ImGui_ImplGlfw_Shutdown();
        ImPlot::DestroyContext();
        ImGui::DestroyContext();
        gui_started = false;
        glfwDestroyWindow(window);
        window = nullptr;
        glfwTerminate();
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "AstraFlow GUI: " << e.what() << '\n';
        if (gui_started) {
            ImGui_ImplOpenGL3_Shutdown();
            ImGui_ImplGlfw_Shutdown();
            ImPlot::DestroyContext();
            ImGui::DestroyContext();
        }
        if (window)
            glfwDestroyWindow(window);
        glfwTerminate();
        return 1;
    }
}

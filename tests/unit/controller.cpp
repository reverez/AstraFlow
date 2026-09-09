#include "astraflow/core/controller.hpp"
#include "astraflow/analysis/fields.hpp"
#include <catch2/catch_test_macros.hpp>
using namespace astraflow;
TEST_CASE("Worker run pause step reset and invalid regeneration", "[unit]") {
    Config c;
    c.backend = "cpu";
    c.precision = "double";
    c.geometry.nx = 16;
    c.geometry.nr = 8;
    c.max_iterations = 100000;
    c.residual_tolerance = 0;
    Controller controller(c);
    auto timeout = std::chrono::milliseconds(5000);
    REQUIRE(controller.wait_for(RunState::Ready, timeout));
    controller.single_step();
    REQUIRE(controller.wait_iterations(1, timeout));
    REQUIRE(controller.wait_for(RunState::Paused, timeout));
    controller.run();
    REQUIRE(controller.wait_iterations(10, timeout));
    controller.pause();
    REQUIRE(controller.wait_for(RunState::Paused, timeout));
    auto old = controller.snapshot();
    int count = old->stats.iterations;
    controller.single_step();
    REQUIRE(controller.wait_iterations(count + 1, timeout));
    REQUIRE(controller.wait_for(RunState::Paused, timeout));
    REQUIRE(controller.snapshot()->stats.iterations == count + 1);
    auto invalid = c;
    invalid.geometry.throat_radius = -1;
    REQUIRE_THROWS(controller.reset(invalid));
    REQUIRE(controller.state() == RunState::Paused);
    c.geometry.nx = 20;
    controller.reset(c);
    REQUIRE(controller.wait_for(RunState::Ready, timeout));
    REQUIRE(controller.snapshot()->mesh->nx == 20);
    REQUIRE(controller.snapshot()->stats.iterations == 0);
    REQUIRE(old->mesh->nx == 16);
    auto snap = controller.snapshot();
    for (int f = 0; f < 9; ++f) {
        auto values = scalar_field(*snap->mesh, snap->conservative, c.settings.gas, Field(f));
        REQUIRE(values.size() == snap->mesh->cells.size());
        for (double v : values)
            REQUIRE(std::isfinite(v));
    }
}

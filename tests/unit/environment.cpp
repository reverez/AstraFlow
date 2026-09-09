#include <catch2/catch_test_macros.hpp>
#include <limits>
TEST_CASE("Required floating-point representations are available", "[unit]") {
    REQUIRE(std::numeric_limits<float>::is_iec559);
    REQUIRE(std::numeric_limits<double>::is_iec559);
    REQUIRE(sizeof(float) == 4);
    REQUIRE(sizeof(double) == 8);
}

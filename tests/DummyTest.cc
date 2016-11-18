#include <catch.hpp>

namespace gscf {
namespace tests {

TEST_CASE("Dummy test", "[dummy]") {
  SECTION("Dummy") { REQUIRE(true); }
}  // TEST_CASE

}  // namespace test
}  // namespace gscf

// Setup the krims exception system for the tests.
#define KRIMS_INIT_EXCEPTION_SYSTEM
#include <krims/ExceptionSystem.hh>

#define CATCH_CONFIG_RUNNER
#include <catch.hpp>
#include <krims/NumComp.hh>

int main(int argc, char* const argv[]) {
  // Throw in case a numerical comparison fails with very
  // detailed information
  krims::NumCompConstants::default_failure_action =
        krims::NumCompActionType::ThrowVerbose;

  // Run catch:
  int result = Catch::Session().run(argc, argv);
  return result;
}

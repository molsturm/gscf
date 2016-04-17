#include <gscf/version.hh>
#include <iostream>
#include <linalgwrap/version.hh>

int main() {
    std::cout << "gscf version: " << gscf::version::version_string()
              << std::endl
              << "linalgwrap version: " << linalgwrap::version::version_string()
              << std::endl;
}

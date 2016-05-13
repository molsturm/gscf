# Finds and sets up the catch target Upstream::catch
#
# such that just linking against them as a dependency does everything
# automatically.
#

add_library(Upstream::catch INTERFACE IMPORTED)
find_path(catch_INCLUDE_DIR catch.hpp
	PATHS
	$ENV{catch_INCLUDE_DIR}
	${PROJECT_SOURCE_DIR}/../linalgwrap/external/catch/include
	DOC "catch include directory"
)

if ("${catch_INCLUDE_DIR}" STREQUAL "catch_INCLUDE_DIR-NOTFOUND")
	message(FATAL_ERROR "Could not find catch include directory. 
Either disable testing of sturmint by setting STURMINT_ENABLE_TESTS to OFF \
or provide a hint where the catch include file can be found via \
the environment variable catch_INCLUDE_DIR.")
endif()

message(STATUS "Found catch at ${catch_INCLUDE_DIR}/catch.hpp")
include_directories(${catch_INCLUDE_DIR})

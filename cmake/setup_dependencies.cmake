# sets these things
#
# 	GSCF_DEPENDENCIES			everyone needs these libraries
# 	GSCF_DEPENDENCIES_DEBUG		debug mode needs these extras
# 	GSCF_DEPENDENCIES_RELEASE		release mode needs these extras
# 	GSCF_DEPENDENCIES_TEST		tests need these extra libraries
#
#       GSCF_DEFINITIONS			definitions for all compilation
#       GSCF_DEFINITIONS_DEBUG		definitions for debug mode
#       GSCF_DEFINITIONS_RELEASE		definitions for release mode
#       

####################
#-- Empty it all --#
####################
set(GSCF_DEPENDENCIES "")
set(GSCF_DEPENDENCIES_DEBUG "")
set(GSCF_DEPENDENCIES_RELEASE "")
set(GSCF_DEPENDENCIES_TEST "")
set(GSCF_DEFINITIONS "")
set(GSCF_DEFINITIONS_DEBUG "")
set(GSCF_DEFINITIONS_RELEASE "")


##################
#-- linalgwrap --#
##################
include(cmake/findLinalgwrap.cmake)

# Check that all required targets are available.
foreach(build ${DRB_BUILD_TYPES})
	if(NOT TARGET "${linalgwrap_${build}_TARGET}")
		message(FATAL_ERROR "We could not find a ${build} version of linalwrap at this location. \
Either disable building a ${build} version of ${CMAKE_PROJECT_NAME} or else \
rebuild linalgwrap with a ${build} version as well.")
	endif()

	# Add dependencies to appropriate versions of gscf
	set(GSCF_DEPENDENCIES_${build} ${GSCF_DEPENDENCIES_${build}} ${linalgwrap_${build}_TARGET})
endforeach()


##############
#-- catch  --#
##############
if(GSCF_ENABLE_TESTS)
	add_library(catch INTERFACE)
	find_path(catch_INCLUDE_DIR catch.hpp
		PATHS
		$ENV{catch_INCLUDE_DIR}
		${CMAKE_SOURCE_DIR}/../linalgwrap/external/catch/include
		DOC "catch include directory"
	)

	if ("${catch_INCLUDE_DIR}" STREQUAL "catch_INCLUDE_DIR-NOTFOUND")
		message(FATAL_ERROR "Could not find catch include directory. 
Either disable testing of gscf by setting GSCF_ENABLE_TESTS to OFF \
or provide a hint where the catch include file can be found via \
the environment variable catch_INCLUDE_DIR.")
	endif()

	message(STATUS "Found catch at ${catch_INCLUDE_DIR}/catch.hpp")

	target_include_directories(catch INTERFACE ${catch_INCLUDE_DIR})
	set(GSCF_DEPENDENCIES_TEST ${GSCF_DEPENDENCIES_TEST} catch)
endif()

##################
#-- rapidcheck --#
##################
if(GSCF_ENABLE_TESTS)
	find_package(rapidcheck REQUIRED CONFIG 
		PATHS 
		${CMAKE_SOURCE_DIR}/../linalgwrap/build
	)

	message(STATUS "Found rapidcheck config at ${rapidcheck_CONFIG}")

	set(GSCF_DEPENDENCIES_TEST ${GSCF_DEPENDENCIES_TEST} Upstream::rapidcheck)
endif()

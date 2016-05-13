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
if (TARGET ${linalgwrap_DEBUG_TARGET} OR TARGET ${linalgwrap_RELEASE_TARGET})
	# If the targets are already defined elsewhere, we are done:
	message(STATUS "Using linalgwrap library provided by build environment.")
else()
	include(cmake/findLinalgwrap.cmake)
endif()

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
	if (TARGET common_catch)
		MESSAGE(STATUS "Using catch provided by build enviroment.")
		set(GSCF_DEPENDENCIES_TEST ${GSCF_DEPENDENCIES_TEST} common_catch)
	else()
		include(cmake/findCatch.cmake)
		set(GSCF_DEPENDENCIES_TEST ${GSCF_DEPENDENCIES_TEST} Upstream::catch)
	endif()

endif()

##################
#-- rapidcheck --#
##################
if(GSCF_ENABLE_TESTS)
	if (TARGET common_rapidcheck)
		MESSAGE(STATUS "Using rapidcheck provided by build environment.")
		set(GSCF_DEPENDENCIES_TEST ${GSCF_DEPENDENCIES_TEST} rapidcheck)
	else()
		find_package(rapidcheck REQUIRED CONFIG 
			PATHS 
			${gscf_SOURCE_DIR}/../linalgwrap/build
		)

		message(STATUS "Found rapidcheck config at ${rapidcheck_CONFIG}")
		set(GSCF_DEPENDENCIES_TEST ${GSCF_DEPENDENCIES_TEST} Upstream::rapidcheck)
	endif()
endif()

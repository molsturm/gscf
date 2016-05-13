# Finds and sets up the linalgwrap targets
#
#   Upstream::linalgwrap.g      (Debug target)
#   Upstream::linalgwrap        (Release target)
#
# Such that just linking against them as a dependency does everything
# automaatically.
#
# If linalgwrap was not compiled as Debug and Release it may happen that
# one of the above targets is not available.
#
# For convenience and future portability the variables 
#      linalgwrap_DEBUG_TARGET      (Name of Debug target)
#      linalgwrap_RELEASE_TARGET    (Name of Release target)
# are available, but note that their presence does not indicate that
# the targets are actually available.
#

# We need linalgwrap version 0.1.0
find_package(linalgwrap 0.1.0 REQUIRED CONFIG 
	PATHS 
	${PROJECT_SOURCE_DIR}/../linalgwrap/build
)
# Now we found the library. Most of the times that's it and we are done.
# But if we got the linalgwrap from a build directory, then it is very
# likely that the header includes cannot be found like this.

# Find dir containing linalgwrap config:
get_filename_component(linalgwrap_config_dir "${linalgwrap_CONFIG}"  DIRECTORY)

set(linalgwrap_DEBUG_TARGET   "Upstream::linalgwrap.g" 
	CACHE INTERNAL "Target name of debug version of linalgwrap")
set(linalgwrap_RELEASE_TARGET "Upstream::linalgwrap"
	CACHE INTERNAL "Target name of release version of linalgwrap")

# So we try to correct this for all existing targets:
foreach(target ${linalgwrap_DEBUG_TARGET} ${linalgwrap_RELEASE_TARGET})
	if (NOT TARGET ${target})
		continue()
	endif()

	# Check that the includes are present:
	get_target_property(LINALGWRAP_INTERFACE_INCLUDES 
		${target} INTERFACE_INCLUDE_DIRECTORIES)

	# If yes continue
	if(NOT "${LINALGWRAP_INTERFACE_INCLUDES}" 
			STREQUAL "LINALGWRAP_INTERFACE_INCLUDES-NOTFOUND")
		continue()
	endif()

	set(LINALGWRAP_INTERFACE_INCLUDES "")

	# Try to find the include dirctory:
	find_path(linalgwrap_INCLUDE_DIR "linalgwrap/version.hh"
		HINTS
		# If we found a build directory, then this is the 
		# path to the include directory
		${linalgwrap_config_dir}/../src/
		PATHS
		$ENV{linalgwrap_INCLUDE_DIR}
		${PROJECT_SOURCE_DIR}/../linalgwrap/src
		DOC "linalgwrap header include directory"
	)

	# Check that the include directory was found
	if ("${linalgwrap_INCLUDE_DIR}" STREQUAL "linalgwrap_INCLUDE_DIR-NOTFOUND")
		message(FATAL_ERROR "Could not find linalgwrap include directory. 
Please provide a hint using the environment variable linalgwrap_INCLUDE_DIR")
	endif()

	# Append to interface includes:
	set(LINALGWRAP_INTERFACE_INCLUDES ${LINALGWRAP_INTERFACE_INCLUDES} ${linalgwrap_INCLUDE_DIR})

	# Check that the linalgwrap/version_defs.hh file can be found in this include
	# directory.
	if(NOT EXISTS "${linalgwrap_INCLUDE_DIR}/linalgwrap/version_defs.hh")
		if(EXISTS "${linalgwrap_config_dir}/src/linalgwrap/version_defs.hh")
			set(LINALGWRAP_INTERFACE_INCLUDES ${LINALGWRAP_INTERFACE_INCLUDES} 
				"${linalgwrap_config_dir}/src"
			)
		else()
			message(FATAL_ERROR "Could not find linalgwrap version_defs.hh file")
		endif()
	endif()

	# Set the interface includes:
	set_target_properties(${target}
		PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${LINALGWRAP_INTERFACE_INCLUDES}"
	)
endforeach()

message(STATUS "Found linalgwrap config at ${linalgwrap_CONFIG}")

## ---------------------------------------------------------------------
## Put licence text here
## ---------------------------------------------------------------------

# Finds and sets up linalgwrap under the target names stored in
#      linalgwrap_DEBUG_TARGET     (Debug version)
#      linalgwrap_RELEASE_TARGET   (Release version)
# such that just linking against it as a dependency does everything
# automatically.
#
# In case the linalgwrap library is not found and AUTOCHECKOUT_MISSING_LIBS is set to
# on, linalgwrap is automatically checked out and built.
# Otherwise a fatal error is produced.
#

#
# Options and properties required
#
option(AUTOCHECKOUT_MISSING_REPOS "Automatically checkout missing repositories" OFF)

#
# -------
#

if (TARGET "${linalgwrap_DEBUG_TARGET}"  OR TARGET "${linalgwrap_RELEASE_TARGET}")
	message(STATUS "Found linalgwrap targets, assume linalgwrap already configured for build.")
	return()
endif()

# Try to find linalgwrap somewhere
find_package(linalgwrap ${LINALGWRAP_VERSION} QUIET CONFIG)
string(TOUPPER "${PROJECT_NAME}" PROJECT_UPPER)
if ("${linalgwrap_DIR}" STREQUAL "linalgwrap_DIR-NOTFOUND")
	if (AUTOCHECKOUT_MISSING_REPOS)
		execute_process(
			COMMAND "sh" "get_linalgwrap.sh"
			WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/external"
			RESULT_VARIABLE RES
		)
		if (NOT RES EQUAL 0)
			message(FATAL_ERROR "Getting linalgwrap from git failed with error: ${RES}")
		endif()

		#
		# Proceed to configure linalgwrap
		#
		add_subdirectory(${PROJECT_SOURCE_DIR}/external/linalgwrap)
		include_directories(${PROJECT_SOURCE_DIR}/external/linalgwrap/src)

		# Extract version from CMakeLists.txt:
		file(STRINGS "${PROJECT_SOURCE_DIR}/external/linalgwrap/CMakeLists.txt"
			VERSION_RAW
			REGEX "linalgwrap VERSION [0-9.]+"
			LIMIT_COUNT 1)
		string(REGEX MATCH "[0-9.]+" GOT_VERSION "${VERSION_RAW}")

		# Compare against what is needed
		if("${GOT_VERSION}" VERSION_LESS "${LINALGWRAP_VERSION}")
			message(FATAL_ERROR "Inconsistency in the repo: \
Version ${LINALGWRAP_VERSION} of linalgwrap was requested, but only version ${GOT_VERSION} \
was found.")
		endif()

		return()
	endif()

	message(FATAL_ERROR "Could not find linalgwrap library.
Either provide the installation prefix of the linalgwrap library in the environment \
variable linalgwrap_DIR or enable autocheckout via -DAUTOCHECKOUT_MISSING_REPOS=ON.")
endif()

#TODO test and remove (same below in #old#)
message(WARNING "This part of findLinalgwrap has never been tested.")

# Setup library targets
set(linalgwrap_DEBUG_TARGET   "Upstream::linalgwrap.g"
	CACHE INTERNAL "Target name of debug version of linalgwrap")
set(linalgwrap_RELEASE_TARGET "Upstream::linalgwrap"
	CACHE INTERNAL "Target name of release version of linalgwrap")

# Check that all required targets are available.
foreach(build ${DRB_BUILD_TYPES})
	if(NOT TARGET "${linalgwrap_${build}_TARGET}")
		message(FATAL_ERROR "We could not find a ${build} version of linalgwrap at this location. \
		Either disable building a ${build} version of ${CMAKE_PROJECT_NAME} or else \
		rebuild linalgwrap with a ${build} version as well.")
	endif()
endforeach()

#TODO check that we don't need this extra stuff:
#old## Now we found the library. Most of the times that's it and we are done.
#old## But if we got the linalgwrap from a build directory, then it is very
#old## likely that the header includes cannot be found like this.
#old#
#old## Find dir containing linalgwrap config:
#old#get_filename_component(linalgwrap_config_dir "${linalgwrap_CONFIG}"  DIRECTORY)
#old#
#old#set(linalgwrap_DEBUG_TARGET   "Upstream::linalgwrap.g" 
#old#	CACHE INTERNAL "Target name of debug version of linalgwrap")
#old#set(linalgwrap_RELEASE_TARGET "Upstream::linalgwrap"
#old#	CACHE INTERNAL "Target name of release version of linalgwrap")
#old#
#old## So we try to correct this for all existing targets:
#old#foreach(target ${linalgwrap_DEBUG_TARGET} ${linalgwrap_RELEASE_TARGET})
#old#	if (NOT TARGET ${target})
#old#		continue()
#old#	endif()
#old#
#old#	# Check that the includes are present:
#old#	get_target_property(LINALGWRAP_INTERFACE_INCLUDES 
#old#		${target} INTERFACE_INCLUDE_DIRECTORIES)
#old#
#old#	# If yes continue
#old#	if(NOT "${LINALGWRAP_INTERFACE_INCLUDES}" 
#old#			STREQUAL "LINALGWRAP_INTERFACE_INCLUDES-NOTFOUND")
#old#		continue()
#old#	endif()
#old#
#old#	set(LINALGWRAP_INTERFACE_INCLUDES "")
#old#
#old#	# Try to find the include dirctory:
#old#	find_path(linalgwrap_INCLUDE_DIR "linalgwrap/version.hh"
#old#		HINTS
#old#		# If we found a build directory, then this is the 
#old#		# path to the include directory
#old#		${linalgwrap_config_dir}/../src/
#old#		PATHS
#old#		$ENV{linalgwrap_INCLUDE_DIR}
#old#		${PROJECT_SOURCE_DIR}/../linalgwrap/src
#old#		DOC "linalgwrap header include directory"
#old#	)
#old#
#old#	# Check that the include directory was found
#old#	if ("${linalgwrap_INCLUDE_DIR}" STREQUAL "linalgwrap_INCLUDE_DIR-NOTFOUND")
#old#		message(FATAL_ERROR "Could not find linalgwrap include directory. 
#old#Please provide a hint using the environment variable linalgwrap_INCLUDE_DIR")
#old#	endif()
#old#
#old#	# Append to interface includes:
#old#	set(LINALGWRAP_INTERFACE_INCLUDES ${LINALGWRAP_INTERFACE_INCLUDES} ${linalgwrap_INCLUDE_DIR})
#old#
#old#	# Check that the linalgwrap/version_defs.hh file can be found in this include
#old#	# directory.
#old#	if(NOT EXISTS "${linalgwrap_INCLUDE_DIR}/linalgwrap/version_defs.hh")
#old#		if(EXISTS "${linalgwrap_config_dir}/src/linalgwrap/version_defs.hh")
#old#			set(LINALGWRAP_INTERFACE_INCLUDES ${LINALGWRAP_INTERFACE_INCLUDES} 
#old#				"${linalgwrap_config_dir}/src"
#old#			)
#old#		else()
#old#			message(FATAL_ERROR "Could not find linalgwrap version_defs.hh file")
#old#		endif()
#old#	endif()
#old#
#old#	# Set the interface includes:
#old#	set_target_properties(${target}
#old#		PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${LINALGWRAP_INTERFACE_INCLUDES}"
#old#	)
#old#endforeach()

message(STATUS "Found linalgwrap config at ${linalgwrap_CONFIG}")

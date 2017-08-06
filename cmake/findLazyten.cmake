## ---------------------------------------------------------------------
##
## Copyright (C) 2017 by the gscf authors
##
## This file is part of gscf.
##
## gscf is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published
## by the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## gscf is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with gscf. If not, see <http://www.gnu.org/licenses/>.
##
## ---------------------------------------------------------------------

# Finds and sets up lazyten under the target names stored in
#      lazyten_DEBUG_TARGET     (Debug version)
#      lazyten_RELEASE_TARGET   (Release version)
# such that just linking against it as a dependency does everything
# automatically.
#
# In case the lazyten library is not found and AUTOCHECKOUT_MISSING_LIBS is set to
# on, lazyten is automatically checked out and built.
# Otherwise a fatal error is produced.
#

#
# Options and properties required
#
option(AUTOCHECKOUT_MISSING_REPOS "Automatically checkout missing repositories" OFF)

#
# -------
#

if (TARGET "${lazyten_DEBUG_TARGET}"  OR TARGET "${lazyten_RELEASE_TARGET}")
	message(STATUS "Found lazyten targets, assume lazyten already configured for build.")
	return()
endif()

# Try to find lazyten somewhere
find_package(lazyten ${LAZYTEN_VERSION} QUIET CONFIG)
mark_as_advanced(lazyten_DIR)

if ("${lazyten_DIR}" STREQUAL "lazyten_DIR-NOTFOUND")
	if (AUTOCHECKOUT_MISSING_REPOS)
		execute_process(
			COMMAND "sh" "get_lazyten.sh"
			WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/external"
			RESULT_VARIABLE RES
		)
		if (NOT RES EQUAL 0)
			message(FATAL_ERROR "Getting lazyten from git failed with error: ${RES}")
		endif()

		#
		# Proceed to configure lazyten
		#
		add_subdirectory(${PROJECT_SOURCE_DIR}/external/lazyten)
		include_directories(${PROJECT_SOURCE_DIR}/external/lazyten/src)
		include_directories(${PROJECT_BINARY_DIR}/external/lazyten/src)

		# Extract version from CMakeLists.txt:
		file(STRINGS "${PROJECT_SOURCE_DIR}/external/lazyten/CMakeLists.txt"
			VERSION_RAW
			REGEX "lazyten VERSION [0-9.]+"
			LIMIT_COUNT 1)
		string(REGEX MATCH "[0-9.]+" GOT_VERSION "${VERSION_RAW}")

		# Compare against what is needed
		if("${GOT_VERSION}" VERSION_LESS "${LAZYTEN_VERSION}")
			message(FATAL_ERROR "Inconsistency in the repo: \
Version ${LAZYTEN_VERSION} of lazyten was requested, but only version ${GOT_VERSION} \
was found.")
		endif()

		return()
	endif()

	message(FATAL_ERROR "Could not find lazyten library.
Either provide the installation prefix of the lazyten library in the environment \
variable lazyten_DIR or enable autocheckout via -DAUTOCHECKOUT_MISSING_REPOS=ON.")
endif()

#TODO test and remove (same below in #old#)
message(WARNING "This part of findLazyten has never been tested.")

# Setup library targets
set(lazyten_DEBUG_TARGET   "Upstream::lazyten.g"
	CACHE INTERNAL "Target name of debug version of lazyten")
set(lazyten_RELEASE_TARGET "Upstream::lazyten"
	CACHE INTERNAL "Target name of release version of lazyten")

# Check that all required targets are available.
foreach(build ${DRB_BUILD_TYPES})
	if(NOT TARGET "${lazyten_${build}_TARGET}")
		message(FATAL_ERROR "We could not find a ${build} version of lazyten at this location. \
		Either disable building a ${build} version of ${CMAKE_PROJECT_NAME} or else \
		rebuild lazyten with a ${build} version as well.")
	endif()
endforeach()

#TODO check that we don't need this extra stuff:
#old## Now we found the library. Most of the times that's it and we are done.
#old## But if we got the lazyten from a build directory, then it is very
#old## likely that the header includes cannot be found like this.
#old#
#old## Find dir containing lazyten config:
#old#get_filename_component(lazyten_config_dir "${lazyten_CONFIG}"  DIRECTORY)
#old#
#old#set(lazyten_DEBUG_TARGET   "Upstream::lazyten.g" 
#old#	CACHE INTERNAL "Target name of debug version of lazyten")
#old#set(lazyten_RELEASE_TARGET "Upstream::lazyten"
#old#	CACHE INTERNAL "Target name of release version of lazyten")
#old#
#old## So we try to correct this for all existing targets:
#old#foreach(target ${lazyten_DEBUG_TARGET} ${lazyten_RELEASE_TARGET})
#old#	if (NOT TARGET ${target})
#old#		continue()
#old#	endif()
#old#
#old#	# Check that the includes are present:
#old#	get_target_property(LAZYTEN_INTERFACE_INCLUDES 
#old#		${target} INTERFACE_INCLUDE_DIRECTORIES)
#old#
#old#	# If yes continue
#old#	if(NOT "${LAZYTEN_INTERFACE_INCLUDES}" 
#old#			STREQUAL "LAZYTEN_INTERFACE_INCLUDES-NOTFOUND")
#old#		continue()
#old#	endif()
#old#
#old#	set(LAZYTEN_INTERFACE_INCLUDES "")
#old#
#old#	# Try to find the include dirctory:
#old#	find_path(lazyten_INCLUDE_DIR "lazyten/version.hh"
#old#		HINTS
#old#		# If we found a build directory, then this is the 
#old#		# path to the include directory
#old#		${lazyten_config_dir}/../src/
#old#		PATHS
#old#		$ENV{lazyten_INCLUDE_DIR}
#old#		${PROJECT_SOURCE_DIR}/../lazyten/src
#old#		DOC "lazyten header include directory"
#old#	)
#old#
#old#	# Check that the include directory was found
#old#	if ("${lazyten_INCLUDE_DIR}" STREQUAL "lazyten_INCLUDE_DIR-NOTFOUND")
#old#		message(FATAL_ERROR "Could not find lazyten include directory. 
#old#Please provide a hint using the environment variable lazyten_INCLUDE_DIR")
#old#	endif()
#old#
#old#	# Append to interface includes:
#old#	set(LAZYTEN_INTERFACE_INCLUDES ${LAZYTEN_INTERFACE_INCLUDES} ${lazyten_INCLUDE_DIR})
#old#
#old#	# Check that the lazyten/version_defs.hh file can be found in this include
#old#	# directory.
#old#	if(NOT EXISTS "${lazyten_INCLUDE_DIR}/lazyten/version_defs.hh")
#old#		if(EXISTS "${lazyten_config_dir}/src/lazyten/version_defs.hh")
#old#			set(LAZYTEN_INTERFACE_INCLUDES ${LAZYTEN_INTERFACE_INCLUDES} 
#old#				"${lazyten_config_dir}/src"
#old#			)
#old#		else()
#old#			message(FATAL_ERROR "Could not find lazyten version_defs.hh file")
#old#		endif()
#old#	endif()
#old#
#old#	# Set the interface includes:
#old#	set_target_properties(${target}
#old#		PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${LAZYTEN_INTERFACE_INCLUDES}"
#old#	)
#old#endforeach()

message(STATUS "Found lazyten config at ${lazyten_CONFIG}")

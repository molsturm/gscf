## ---------------------------------------------------------------------
##
## Copyright (C) 2016 by the gscf authors
##
## This file is part of gscf.
##
## gscf is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published
## by the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## gscf is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU Lesser General Public License for more details.
##
## You should have received a copy of the GNU Lesser General Public License
## along with gscf. If not, see <http://www.gnu.org/licenses/>.
##
## ---------------------------------------------------------------------

cmake_minimum_required(VERSION 3.0.0)

# Load the DebugReleaseBuild module for a simultaneous Debug and Release
# and release build and initialise it.
include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/IncludeKrimsCmakeModule.cmake")
include_krims_cmake_module(DebugReleaseBuild)
drb_init()

##########################################################################
# Global options

# Path where cmake wants the modules and config packages:
set(PackageModuleLocation share/cmake)

# setup global options (needs to be done before call to project())
option(GSCF_ENABLE_TESTS "Build gscf tests" ON)
option(GSCF_ENABLE_EXAMPLES "Build gscf examples" ON)

##########################################################################
# Project setup

# enable the new version feature introduced in cmake 3.0.0
cmake_policy(SET CMP0048 NEW)
project(gscf VERSION 0.0.0)
message("\n# Configuring ${PROJECT_NAME}\n#")

# Global include directories for this project:
include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories("${PROJECT_BINARY_DIR}/src")

# enable testing of this project
enable_testing()

##########################################################################
# Setup hard and optional dependencies and find components

# Determine and setup compiler flags, enforcing C++11
drb_setup_compiler_flags(11)

# setup hard dependencies
include(cmake/setup_dependencies.cmake)

# setup optional dependencies and features.
include(cmake/setup_optional.cmake)

##########################################################################
# Pass on to subdirectories:

# The src directory does most of the work: Here the library is build and
# installation is configured.
add_subdirectory(src)

# Add subdirectories for the test.
if(GSCF_ENABLE_TESTS)
	add_subdirectory(tests)
endif()

# Add subdirectories for the examples.
if(GSCF_ENABLE_EXAMPLES)
	add_subdirectory(examples)
endif()

##########################################################################
# Finally install package:
include(cmake/install_package.cmake)

# and print a closing message
message("#\n# Configuring ${PROJECT_NAME} done")

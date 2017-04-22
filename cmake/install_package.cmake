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

# Installs the cmake apckage information this project provides
#
# Requires the variable PackageModuleLocation to be set.

# Write a basic version file for gscf
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
	"${gscf_BINARY_DIR}/gscfConfigVersion.cmake"
	COMPATIBILITY AnyNewerVersion
)

# Adjust a configure file
configure_file(cmake/gscfConfig.cmake.in
	"${gscf_BINARY_DIR}/gscfConfig.cmake"
	COPYONLY
)

# Set an export location:
install(FILES
	"${gscf_BINARY_DIR}/gscfConfig.cmake"
	"${gscf_BINARY_DIR}/gscfConfigVersion.cmake"
	DESTINATION "${PackageModuleLocation}/gscf"
	COMPONENT devel
)


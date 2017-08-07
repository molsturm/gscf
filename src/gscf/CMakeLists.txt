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

#
# set the source files.
#
set(GSCF_SOURCES
	ScfBaseKeys.cc
	PulayDiisScfKeys.cc
	TruncatedOptDampScfKeys.cc
	version.cc
)

# Dump the current version as an include file.
# Write the config file config.hh (into the binary dir)
configure_file("config.hh.in"  "config.hh")

# Add the libraries for Debug and Release builds:
drb_add_library(gscf
	DBGSUFFIX  ".g"
	FILES ${GSCF_SOURCES}
)

drb_set_target_properties(
	ALL gscf
	PROPERTIES
	VERSION "${PROJECT_VERSION}"
)

drb_target_link_libraries(ALL     gscf ${GSCF_DEPENDENCIES})
drb_target_link_libraries(DEBUG   gscf ${GSCF_DEPENDENCIES_DEBUG})
drb_target_link_libraries(RELEASE gscf ${GSCF_DEPENDENCIES_RELEASE})

# set the installation property for these targets:
install(TARGETS ${gscf_TARGETS}
	EXPORT GscfTargets
	DESTINATION lib
	COMPONENT library
)

# install the includes:
install(DIRECTORY .
	DESTINATION "include/gscf"
	COMPONENT devel
	FILES_MATCHING PATTERN "*.hh"
)
install(FILES "${gscf_BINARY_DIR}/src/gscf/config.hh"
	DESTINATION "include/gscf"
	COMPONENT devel
)

# Export the target specifications for gscf.
export(EXPORT GscfTargets
	FILE "${gscf_BINARY_DIR}/gscfTargets.cmake"
	NAMESPACE Upstream::
)

# and install them appropriately:
install(EXPORT GscfTargets
	FILE "gscfTargets.cmake"
	NAMESPACE Upstream::
	DESTINATION "${PackageModuleLocation}/gscf"
)


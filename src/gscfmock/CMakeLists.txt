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

# The libraries for Debug and Release builds:
drb_add_library(gscfmock
	DBGSUFFIX  ".g"
	FILES
	Integrals/IntegralsSturmian14.cc
	FockMatrix.cc
	pulay_error.cc
)

drb_set_target_properties(
	ALL gscfmock
	PROPERTIES
	VERSION "${PROJECT_VERSION}"
)

drb_target_compile_definitions(ALL     gscfmock PUBLIC ${GSCF_DEFINITIONS})
drb_target_compile_definitions(DEBUG   gscfmock PUBLIC ${GSCF_DEFINITIONS_DEBUG})
drb_target_compile_definitions(RELEASE gscfmock PUBLIC ${GSCF_DEFINITIONS_RELEASE})

drb_target_link_libraries(ALL     gscfmock ${GSCF_DEPENDENCIES})
drb_target_link_libraries(DEBUG   gscfmock ${GSCF_DEPENDENCIES_DEBUG})
drb_target_link_libraries(RELEASE gscfmock ${GSCF_DEPENDENCIES_RELEASE})

# Note:
# Since these files are only needed for the examples
# and the tests we do not install them at all.

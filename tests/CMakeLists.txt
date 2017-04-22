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
# The sources for the test executable
#
set(GSCF_TESTS_SOURCES
	FunctionalityTest.cc
	main.cc
)

drb_add_executable(gscf_tests
	DBGSUFFIX  "_DEBUG"
	RELSUFFIX  "_RELEASE" 
	FILES ${GSCF_TESTS_SOURCES}
)

drb_target_compile_definitions(ALL     gscf_tests 
	PRIVATE ${GSCF_DEFINITIONS})
drb_target_compile_definitions(DEBUG   gscf_tests 
	PRIVATE ${GSCF_DEFINITIONS_DEBUG})
drb_target_compile_definitions(RELEASE gscf_tests 
	PRIVATE ${GSCF_DEFINITIONS_RELEASE})

drb_target_link_libraries(ALL     gscf_tests 
	${GSCF_DEPENDENCIES} ${GSCF_DEPENDENCIES_TEST})
drb_target_link_libraries(DEBUG   gscf_tests 
	${gscf_DEBUG_TARGET} ${gscfmock_DEBUG_TARGET} ${GSCF_DEPENDENCIES_DEBUG})
drb_target_link_libraries(RELEASE gscf_tests 
	${gscf_RELEASE_TARGET} ${gscfmock_RELEASE_TARGET} ${GSCF_DEPENDENCIES_RELEASE})

# Add the test executables
FOREACH(target ${gscf_tests_TARGETS})
	add_test(NAME ${target} COMMAND ${target})
ENDFOREACH()


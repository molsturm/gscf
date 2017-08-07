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

# sets these things
#
#       GSCF_DEPENDENCIES			everyone needs these libraries
#       GSCF_DEPENDENCIES_DEBUG			debug mode needs these extras
#       GSCF_DEPENDENCIES_RELEASE		release mode needs these extras
#       GSCF_DEPENDENCIES_TEST			tests need these extra libraries
#

####################
#-- Empty it all --#
####################
set(GSCF_DEPENDENCIES "")
set(GSCF_DEPENDENCIES_DEBUG "")
set(GSCF_DEPENDENCIES_RELEASE "")
set(GSCF_DEPENDENCIES_TEST "")

############################
#-- rapidcheck and catch --#
############################
if(GSCF_ENABLE_TESTS)
	# We need to setup rapidcheck and catch for the tests:
	include(cmake/findRapidcheck.cmake)
	set(GSCF_DEPENDENCIES_TEST ${GSCF_DEPENDENCIES_TEST} ${rapidcheck_TARGET})

	include(cmake/findCatch.cmake)
	set(GSCF_DEPENDENCIES_TEST ${GSCF_DEPENDENCIES_TEST} ${catch_TARGET})
endif()

#############
#-- krims --#
#############
set(KRIMS_VERSION 0.1.0)
include(cmake/findKrims.cmake)

foreach (build ${DRB_BUILD_TYPES})
	set(GSCF_DEPENDENCIES_${build} ${GSCF_DEPENDENCIES_${build}} ${krims_${build}_TARGET})
endforeach()

###############
#-- lazyten --#
###############
set(LAZYTEN_VERSION 0.3.0)
include(cmake/findLazyten.cmake)

foreach (build ${DRB_BUILD_TYPES})
	set(GSCF_DEPENDENCIES_${build} ${GSCF_DEPENDENCIES_${build}} ${lazyten_${build}_TARGET})
endforeach()

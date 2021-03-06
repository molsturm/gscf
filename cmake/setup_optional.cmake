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

# Setup optional dependencies and features
# alters these things
#
#       GSCF_DEPENDENCIES		everyone needs these libraries
#       GSCF_DEPENDENCIES_DEBUG		debug mode needs these extras
#       GSCF_DEPENDENCIES_RELEASE	release mode needs these extras
#       GSCF_DEPENDENCIES_TEST		tests need these extra libraries
#

####################
#-- C++ standard --#
####################
if (NOT CMAKE_CXX_STANDARD VERSION_LESS 14)
	message(STATUS "Detected C++14 support: Setting GSCF_HAVE_CXX14")
	set(GSCF_HAVE_CXX14 ON)
endif()
if (NOT CMAKE_CXX_STANDARD VERSION_LESS 17)
	message(STATUS "Detected C++17 support: Setting GSCF_HAVE_CXX17")
	set(GSCF_HAVE_CXX17 ON)
endif()


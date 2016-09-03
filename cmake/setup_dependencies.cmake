# sets these things
#
# 	GSCF_DEPENDENCIES			everyone needs these libraries
# 	GSCF_DEPENDENCIES_DEBUG		debug mode needs these extras
# 	GSCF_DEPENDENCIES_RELEASE		release mode needs these extras
# 	GSCF_DEPENDENCIES_TEST		tests need these extra libraries
#
#       GSCF_DEFINITIONS			definitions for all compilation
#       GSCF_DEFINITIONS_DEBUG		definitions for debug mode
#       GSCF_DEFINITIONS_RELEASE		definitions for release mode
#       

####################
#-- Empty it all --#
####################
set(GSCF_DEPENDENCIES "")
set(GSCF_DEPENDENCIES_DEBUG "")
set(GSCF_DEPENDENCIES_RELEASE "")
set(GSCF_DEPENDENCIES_TEST "")
set(GSCF_DEFINITIONS "")
set(GSCF_DEFINITIONS_DEBUG "")
set(GSCF_DEFINITIONS_RELEASE "")

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
# Find at least version 0.0.0
set(KRIMS_VERSION 0.0.0)
include(cmake/findKrims.cmake)

foreach (build ${DRB_BUILD_TYPES})
	set(LINALGWRAP_DEPENDENCIES_${build} ${LINALGWRAP_DEPENDENCIES_${build}} ${krims_${build}_TARGET})
endforeach()

##################
#-- linalgwrap --#
##################
# Find at least version 0.2.0
set(LINALGWRAP_VERSION 0.2.0)
include(cmake/findLinalgwrap.cmake)

foreach (build ${DRB_BUILD_TYPES})
	set(GSCF_DEPENDENCIES_${build} ${GSCF_DEPENDENCIES_${build}} ${linalgwrap_${build}_TARGET})
endforeach()

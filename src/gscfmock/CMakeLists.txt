# The libraries for Debug and Release builds:
drb_add_library(gscfmock
	DBGSUFFIX  ".g"
	FILES
	#
	Integrals/IntegralsSturmian14.cc
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

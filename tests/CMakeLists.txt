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


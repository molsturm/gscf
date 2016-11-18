#
# The sources for this example executable
#
# Compile and link:
add_executable(dummy_scf main.cc)

set_target_properties(
	dummy_scf
	PROPERTIES
	LINK_FLAGS "${CMAKE_EXE_LINKER_FLAGS_${EXAMPLES_BUILD_TYPE}}" #${CMAKE_EXE_LINKER_FLAGS} are added
	LINKER_LANGUAGE "CXX"
	COMPILE_FLAGS "${CMAKE_CXX_FLAGS_${EXAMPLES_BUILD_TYPE}}" #${CMAKE_CXX_FLAGS} are added automaticaly
	COMPILE_DEFINITIONS "${EXAMPLES_BUILD_TYPE}"
)

# setup definitions.
target_compile_definitions(dummy_scf PRIVATE "${GSCF_DEFINITIONS};${GSCF_DEFINITIONS_${EXAMPLES_BUILD_TYPE}}")

# link it with the appropriate gscf and gscfmock library target
target_link_libraries(dummy_scf
	${gscf_${EXAMPLES_BUILD_TYPE}_TARGET}
	${gscfmock_${EXAMPLES_BUILD_TYPE}_TARGET}
	${GSCF_DEPENDENCIES}
	${GSCF_DEPENDENCIES_${EXAMPLES_BUILD_TYPE}}
)

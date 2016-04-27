#
# The sources for this example executable
#
set(EXAMPLE_SCF_DUMMY_SOURCES
	integrals/IntegralsSturmian14.cc
	main.cc)

# Compile and link:
add_executable(scf_dummy ${EXAMPLE_SCF_DUMMY_SOURCES})

set_target_properties(
	scf_dummy
	PROPERTIES
	LINK_FLAGS "${CMAKE_EXE_LINKER_FLAGS_${EXAMPLES_BUILD_TYPE}}" #${CMAKE_EXE_LINKER_FLAGS} are added
	LINKER_LANGUAGE "CXX"
	COMPILE_FLAGS "${CMAKE_CXX_FLAGS_${EXAMPLES_BUILD_TYPE}}" #${CMAKE_CXX_FLAGS} are added automaticaly
	COMPILE_DEFINITIONS "${EXAMPLES_BUILD_TYPE}"
)

# setup definitions.
target_compile_definitions(scf_dummy PRIVATE "${GSCF_DEFINITIONS};${GSCF_DEFINITIONS_${EXAMPLES_BUILD_TYPE}}")

# link it with the appropriate linalgwrap library target
target_link_libraries(scf_dummy
	${gscf_${EXAMPLES_BUILD_TYPE}_TARGET}
	${GSCF_DEPENDENCIES}
	${GSCF_DEPENDENCIES_${EXAMPLES_BUILD_TYPE}}
)

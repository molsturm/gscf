#
# The sources for this example executable
#
set(EXAMPLE_ATOMSCF_SOURCES
	main.cc)

# Compile and link:
add_executable(atom_scf ${EXAMPLE_ATOMSCF_SOURCES})

set_target_properties(
	atom_scf
	PROPERTIES
	LINK_FLAGS "${CMAKE_EXE_LINKER_FLAGS_${EXAMPLES_BUILD_TYPE}}" #${CMAKE_EXE_LINKER_FLAGS} are added
	LINKER_LANGUAGE "CXX"
	COMPILE_FLAGS "${CMAKE_CXX_FLAGS_${EXAMPLES_BUILD_TYPE}}" #${CMAKE_CXX_FLAGS} are added automaticaly
)

# setup definitions.
target_compile_definitions(atom_scf PRIVATE "${GSCF_DEFINITIONS};${GSCF_DEFINITIONS_${EXAMPLES_BUILD_TYPE}}")

# link it with the appropriate linalgwrap library target
target_link_libraries(atom_scf
	${gscf_${EXAMPLES_BUILD_TYPE}_TARGET}
	${GSCF_DEPENDENCIES}
	${GSCF_DEPENDENCIES_${EXAMPLES_BUILD_TYPE}}
)

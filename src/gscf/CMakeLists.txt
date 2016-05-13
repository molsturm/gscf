#
# set the source files.
#
set(GSCF_SOURCES
	version.cc
)

# Dump the current version as an include file.
include_linalgwrap_cmake_module(WriteVersionHeader)
write_version_header("${CMAKE_BINARY_DIR}/src/gscf/version_defs.hh")

# Add the libraries for Debug and Release builds:
drb_add_library(gscf
	DBGSUFFIX  ".g"
	FILES ${GSCF_SOURCES}
)

drb_set_target_properties(
	ALL gscf
	PROPERTIES
	VERSION "${PROJECT_VERSION}"
)

drb_target_compile_definitions(ALL     gscf PUBLIC ${GSCF_DEFINITIONS})
drb_target_compile_definitions(DEBUG   gscf PUBLIC ${GSCF_DEFINITIONS_DEBUG})
drb_target_compile_definitions(RELEASE gscf PUBLIC ${GSCF_DEFINITIONS_RELEASE})

drb_target_link_libraries(ALL     gscf ${GSCF_DEPENDENCIES})
drb_target_link_libraries(DEBUG   gscf ${GSCF_DEPENDENCIES_DEBUG})
drb_target_link_libraries(RELEASE gscf ${GSCF_DEPENDENCIES_RELEASE})

# set the installation property for these targets:
install(TARGETS ${gscf_TARGETS}
	EXPORT GscfTargets
	DESTINATION lib 
	COMPONENT library
)

# install the includes:
install(DIRECTORY . 
	DESTINATION "include/gscf"
	COMPONENT devel
	FILES_MATCHING PATTERN "*.hh"
)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/version_defs.hh"
	DESTINATION "include/gscf"
	COMPONENT devel
)

# Export the target specifications for gscf.
export(EXPORT GscfTargets
	FILE "${gscf_BINARY_DIR}/gscfTargets.cmake"
	NAMESPACE Upstream::
)

# and install them appropriately:
install(EXPORT GscfTargets
	FILE "gscfTargets.cmake"
	NAMESPACE Upstream::
	DESTINATION "${PackageModuleLocation}/gscf"
)


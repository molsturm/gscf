# Installs the cmake apckage information this project provides
#
# Requires the variable PackageModuleLocation to be set.

# Write a basic version file for gscf
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
	"${gscf_BINARY_DIR}/gscfConfigVersion.cmake"
	COMPATIBILITY AnyNewerVersion
)

# Adjust a configure file
configure_file(cmake/gscfConfig.cmake.in
	"${gscf_BINARY_DIR}/gscfConfig.cmake"
	COPYONLY
)

# Set an export location:
install(FILES
	"${gscf_BINARY_DIR}/gscfConfig.cmake"
	"${gscf_BINARY_DIR}/gscfConfigVersion.cmake"
	DESTINATION "${PackageModuleLocation}/gscf"
	COMPONENT devel
)


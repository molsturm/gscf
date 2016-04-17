# Try to find a cmake module which usually gets shipped with linalgwrap.

function(include_linalgwrap_cmake_module MODULE)
	# First try to load it plainly as a module:
	include(${MODULE} OPTIONAL RESULT_VARIABLE RES)

	if ("${RES}" STREQUAL "NOTFOUND")
		# We could not "just" find it.
		# Now try some hints:
		set(ModuleHints 
			# if the user specified the location of the module explictly
			"$ENV{${MODULE}_DIR}/${MODULE}.cmake"
			# if the user provided a hint for linalgwrap:
			"$ENV{linalgwrap_DIR}/share/cmake/modules/${MODULE}.cmake"
			# if linalgwrap is in the same top directory
			"${CMAKE_SOURCE_DIR}/../linalgwrap/cmake/modules/${MODULE}.cmake"
		)

		foreach(hfile ${ModuleHints})
			include(${hfile} OPTIONAL RESULT_VARIABLE RES)
			if (NOT "${RES}" STREQUAL "NOTFOUND")
				message(STATUS "Found ${MODULE} file at ${RES}")
				return()
			endif()
		endforeach()

		message(FATAL_ERROR "Could not find the ${MODULE} module.
Try specifying the enviroment variable ${MODULE}_DIR for a hint toward the directory \
containing the ${MODULE} module. 
Alternatively you can also specify in the env. variable linalgwrap_DIR the installation prefix\
which was used when installing the linalgwrap binaries (which contains this module).")
	else()
		message(STATUS "Using system-provided ${MODULE} file.")
	endif()
endfunction(include_linalgwrap_cmake_module)

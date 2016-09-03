if (BUILD_EXTERNAL_RAPIDCHECK)
	# Check whether rapidcheck is there
	if(NOT EXISTS "${CMAKE_CURRENT_LIST_DIR}/rapidcheck/CMakeLists.txt")
		message(FATAL_ERROR "Could not find rapidcheck at expected place. Try a manual git checkout.")
	endif()

	# Set option such that rapidcheck tests are built.
	set(RC_ENABLE_TESTS ON CACHE BOOL "Build RapidCheck tests")

	# Change compiler flags (CMAKE_CXX_FLAGS) to make fresh build config
	set(CMAKE_CXX_FLAGS_STORED_TMP ${CMAKE_CXX_FLAGS})
	set(CMAKE_CXX_FLAGS "")
	#enable_if_compiles(CMAKE_CXX_FLAGS "-Wno-gnu-zero-variadic-macro-arguments")

	# Add the rapidcheck subdirectory and configure its built.
	add_subdirectory(rapidcheck)

	# undo the changes to the compiler flags:
	set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS_STORED_TMP})
	unset(CMAKE_CXX_FLAGS_STORED_TMP)
endif()

if (BUILD_EXTERNAL_KRIMS)
	# Check whether krims is there
	if(NOT EXISTS "${CMAKE_CURRENT_LIST_DIR}/krims/CMakeLists.txt")
		message(FATAL_ERROR "Could not find krims at expected place. Try a manual git checkout.")
	endif()

	add_subdirectory(krims)
endif()

if (BUILD_EXTERNAL_LINALGWRAP)
	# Check whether linalgwrap is there
	if(NOT EXISTS "${CMAKE_CURRENT_LIST_DIR}/linalgwrap/CMakeLists.txt")
		message(FATAL_ERROR "Could not find linalgwrap at expected place. Try a manual git checkout.")
	endif()

	add_subdirectory(linalgwrap)
endif()

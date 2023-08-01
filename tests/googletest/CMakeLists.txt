include(FetchContent)

if(${CMAKE_VERSION} VERSION_GREATER "3.24")
	# Fixes a warning about DOWNLOAD_EXTRACT_TIMESTAMP.
	# See https://cmake.org/cmake/help/v3.24/policy/CMP0135.html#policy:CMP0135
	cmake_policy(SET CMP0135 NEW)
endif()

FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/6f6ab4212aa02cfe02e480711246da4fc17b0761.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)


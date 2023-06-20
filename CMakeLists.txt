cmake_minimum_required(VERSION 3.14)

project(chemmisol VERSION 0.1)

# GoogleTest requires at least C++14
set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS_DEBUG "\
${CMAKE_CXX_FLAGS_DEBUG} -std=c++14 \
-Wall -Wextra -Wpedantic -Werror -Wno-error=pedantic \
-Woverloaded-virtual -Wnon-virtual-dtor \
-fsanitize=address \
")

include(FetchContent)
FetchContent_Declare(
  argparse
  GIT_REPOSITORY https://github.com/p-ranav/argparse.git
  GIT_TAG        557948f1236db9e27089959de837cc23de6c6bbd
)
FetchContent_MakeAvailable(argparse)

add_library(chemmisol-units-lib
	src/units.cpp
	)
target_include_directories(chemmisol-units-lib PUBLIC include)
add_executable(convert
	src/convert.cpp
	)
target_link_libraries(convert
	chemmisol-units-lib
	argparse::argparse
	)

add_library(chemmisol-lib
	src/chemical.cpp
	src/gauss.cpp
	src/experiment.cpp
	)
target_link_libraries(chemmisol-lib
	PUBLIC chemmisol-units-lib)
target_include_directories(chemmisol-lib PUBLIC include)

add_executable(chemmisol
	src/main.cpp)
target_link_libraries(chemmisol chemmisol-lib)

FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        a3580180d16923d6d5f488e20b3814608a892f17
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)


enable_testing()

add_executable(
	convert_test
	tests/units.cpp
	)
target_link_libraries(
  convert_test
  chemmisol-units-lib
  GTest::gtest_main
  GTest::gmock_main
)

add_executable(
  chemmisol-tests
  tests/linear.cpp
  tests/newton.cpp
  tests/chemical.cpp
  tests/gauss.cpp
)
target_link_libraries(
  chemmisol-tests
  chemmisol-lib
  GTest::gtest_main
  GTest::gmock_main
)

include(GoogleTest)
gtest_discover_tests(convert_test)
gtest_discover_tests(chemmisol-tests)

# Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
	file(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/docs)
	set(DOXYGEN_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/docs)
	#set(DOXYGEN_MACRO_EXPANSION "YES")
	set(DOXYGEN_HTML_OUTPUT ".")
	set(DOXYGEN_SORT_MEMBER_DOCS "NO")
	set(DOXYGEN_HAVE_DOT "YES")
	#set(DOXYGEN_TEMPLATE_RELATIONS "YES")
	doxygen_add_docs(doc
		${CMAKE_CURRENT_SOURCE_DIR}/include
		)
endif()
